package drafts;

import java.util.ArrayList;
import java.util.Arrays;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.MonteCarloSimulationModel;
import net.finmath.montecarlo.RandomVariableFromDoubleArray;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAAD;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import net.finmath.montecarlo.conditionalexpectation.RegressionBasisFunctionsProvider;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.ConditionalExpectationEstimator;
import net.finmath.stochastic.RandomVariable;

public class SimpleLowerBoundEstimation extends AbstractSimpleBoundEstimation
		implements RegressionBasisFunctionsProvider {

	public enum BasisFunctionType {
		SwapRates, ForwardRates
	}

	RandomVariable triggerValues;
	RegressionBasisFunctionsProvider regressionBasisFunctionsProvider;

	public SimpleLowerBoundEstimation(BasisFunctionType basisFunctionType) {

		if (basisFunctionType == BasisFunctionType.ForwardRates)
			setRegressionBasisFunctionsProvider(getBasisFunctionsProviderWithSwapRates());
		if (basisFunctionType == BasisFunctionType.SwapRates)
			setRegressionBasisFunctionsProvider(getBasisFunctionsProviderWithForwardRates());

	}

	public SimpleLowerBoundEstimation() {

	}

	@Override
	public RandomVariable getTriggerValues() {
		return triggerValues;
	}

	public void setTriggerValues(RandomVariableDifferentiableAAD triggerValues) {
		this.triggerValues = triggerValues;
	}

	@Override
	protected RandomVariable calculateTriggerValues(int period, double fixingDate,
			LIBORModelMonteCarloSimulationModel model) throws CalculationException {

		SimplestExerciseStrategy exerciseStrategy = new SimplestExerciseStrategy();
		// Calculate the exercise criteria (exercise if the following trigger is
		// negative)
		RandomVariable triggerValuesDiscounted = exerciseStrategy.getExerciseIndicators(continuationValue,
				exerciseValue);

		// Remove foresight through condition expectation
		ConditionalExpectationEstimator conditionalExpectationOperator = getConditionalExpectationEstimator(fixingDate,
				model);

		// Calculate conditional expectation. Note that no discounting (numeraire
		// division) is required!
		triggerValues = triggerValuesDiscounted.getConditionalExpectation(conditionalExpectationOperator);

		// cache values of underlying (including perfect foresight) and E[U_i|F_T_i-1]
		// as needed in 18.24

		cacheValuesOfUnderlying[period] = exerciseValue;
		cacheConditionalExpectations[period] = triggerValuesDiscounted.add(exerciseValue).cache();

		return triggerValues;

	}

	/**
	 * Return the conditional expectation estimator suitable for this product.
	 *
	 * @param fixingDate The condition time.
	 * @param model      The model
	 * @return The conditional expectation estimator suitable for this product
	 * @throws net.finmath.exception.CalculationException Thrown if the valuation
	 *         fails, specific cause may be available via the <code>cause()</code>
	 *         method.
	 */
	public ConditionalExpectationEstimator getConditionalExpectationEstimator(double fixingDate,
			LIBORModelMonteCarloSimulationModel model) throws CalculationException {

		RandomVariable[] regressionBasisFunctions = regressionBasisFunctionsProvider != null
				? regressionBasisFunctionsProvider.getBasisFunctions(fixingDate, model)
				: getBasisFunctions(fixingDate, model);
		MonteCarloConditionalExpectationRegression condExpEstimator = new MonteCarloConditionalExpectationRegression(
				regressionBasisFunctions);

		return condExpEstimator;
	}

	@Override
	public RandomVariable[] getBasisFunctions(double evaluationTime, MonteCarloSimulationModel model)
			throws CalculationException {
		LIBORModelMonteCarloSimulationModel liborModel = (LIBORModelMonteCarloSimulationModel) model;
		return getBasisFunctions(evaluationTime, liborModel);
	}

	/**
	 * Provides a set of \( \mathcal{F}_{t} \)-measurable random variables which can
	 * serve as regression basis functions.
	 *
	 * @param evaluationTime The evaluation time \( t \) at which the basis function
	 *                       should be observed.
	 * @param model          The Monte-Carlo model used to derive the basis
	 *                       function.
	 * @return An \( \mathcal{F}_{t} \)-measurable random variable.
	 * @throws CalculationException Thrown if derivation of the basis function
	 *                              fails.
	 */
	public RandomVariable[] getBasisFunctions(double fixingDate, LIBORModelMonteCarloSimulationModel model)
			throws CalculationException {

		ArrayList<RandomVariable> basisFunctions = new ArrayList<>();
		double[] fixingDates = bermudanOption.getFixingDates();
		double[] paymentDates = bermudanOption.getPaymentDates();
		// Constant
		RandomVariable one = new RandomVariableFromDoubleArray(1.0);
		basisFunctions.add(one);

		int fixingDateIndex = Arrays.binarySearch(fixingDates, fixingDate);
		if (fixingDateIndex < 0) {
			fixingDateIndex = -fixingDateIndex;
		}
		if (fixingDateIndex >= fixingDates.length) {
			fixingDateIndex = fixingDates.length - 1;
		}

		// forward rate to the next period
		RandomVariable rateShort = model.getLIBOR(fixingDate, fixingDate, paymentDates[fixingDateIndex]);
		RandomVariable discountShort = rateShort.mult(paymentDates[fixingDateIndex] - fixingDate).add(1.0).invert();
		basisFunctions.add(discountShort);
		basisFunctions.add(discountShort.pow(2.0));
		// basisFunctions.add(rateShort.pow(3.0));

		// forward rate to the end of the product
		RandomVariable rateLong = model.getLIBOR(fixingDate, fixingDates[fixingDateIndex],
				paymentDates[paymentDates.length - 1]);
		RandomVariable discountLong = rateLong
				.mult(paymentDates[paymentDates.length - 1] - fixingDates[fixingDateIndex]).add(1.0).invert();
		basisFunctions.add(discountLong);
		basisFunctions.add(discountLong.pow(2.0));
		// basisFunctions.add(rateLong.pow(3.0));

		// Numeraire
		RandomVariable numeraire = model.getNumeraire(fixingDate).invert();
		basisFunctions.add(numeraire);
		// basisFunctions.add(numeraire.pow(2.0));
		// basisFunctions.add(numeraire.pow(3.0));

		return basisFunctions.toArray(new RandomVariable[basisFunctions.size()]);
	}

	/*
	 * Some popular variants to create regression basis functions
	 */

	public RegressionBasisFunctionsProvider getBasisFunctionsProviderWithSwapRates() {
		return new RegressionBasisFunctionsProvider() {
			@Override
			public RandomVariableDifferentiableAAD[] getBasisFunctions(double evaluationTime,
					MonteCarloSimulationModel monteCarloModel) throws CalculationException {
				LIBORModelMonteCarloSimulationModel model = (LIBORModelMonteCarloSimulationModel) monteCarloModel;

				double[] exerciseDates = bermudanOption.getFixingDates();
				double[] regressionBasisfunctionTimes = bermudanOption.getPaymentDates();

				ArrayList<RandomVariable> basisFunctions = new ArrayList<>();

				double exerciseTime = evaluationTime;

				int exerciseIndex = Arrays.binarySearch(regressionBasisfunctionTimes, exerciseTime);
				if (exerciseIndex < 0) {
					exerciseIndex = -exerciseIndex;
				}
				if (exerciseIndex >= exerciseDates.length) {
					exerciseIndex = exerciseDates.length - 1;
				}

				// Constant
				RandomVariable one = new RandomVariableFromDoubleArray(1.0);
				RandomVariable basisFunction = one;
				basisFunctions.add(basisFunction);
				// Numeraire (adapted to multicurve framework)
				RandomVariable discountFactor = model.getNumeraire(exerciseTime).invert();
				basisFunctions.add(discountFactor);
				/*
				 * Add swap rates of underlyings.
				 */

				for (int exerciseIndexUnderlying = exerciseIndex; exerciseIndexUnderlying < exerciseDates.length; exerciseIndexUnderlying++) {
					RandomVariable floatLeg = getValueOfLegAnalytic(exerciseTime, exerciseDates,
							regressionBasisfunctionTimes, model, true, 0.0);
					RandomVariable annuity = getValueOfLegAnalytic(exerciseTime, exerciseDates,
							regressionBasisfunctionTimes, model, false, 1.0);
					RandomVariable swapRate = floatLeg.div(annuity);
					basisFunction = swapRate.mult(discountFactor);
					basisFunctions.add(basisFunction);
					basisFunctions.add(basisFunction.squared());
				}

				// forward rate to the next period
				RandomVariable rateShort = model.getLIBOR(exerciseTime, exerciseTime,
						regressionBasisfunctionTimes[exerciseIndex]);
				basisFunctions.add(rateShort);
				basisFunctions.add(rateShort.pow(2.0));

				return basisFunctions.toArray(new RandomVariableDifferentiableAAD[basisFunctions.size()]);
			}
		};
	}

	public RegressionBasisFunctionsProvider getBasisFunctionsProviderWithForwardRates() {
		return new RegressionBasisFunctionsProvider() {
			@Override
			public RandomVariable[] getBasisFunctions(double evaluationTime, MonteCarloSimulationModel monteCarloModel)
					throws CalculationException {
				LIBORModelMonteCarloSimulationModel model = (LIBORModelMonteCarloSimulationModel) monteCarloModel;

				double[] exerciseDates = bermudanOption.getFixingDates();
				double swapEndDate = bermudanOption.getPaymentDates()[exerciseDates.length - 1];
				double[] regressionBasisfunctionTimes = new double[exerciseDates.length + 1];
				int numberOfPeriods = exerciseDates.length;
				for (int i = 0; i < exerciseDates.length; i++)
					regressionBasisfunctionTimes[i] = exerciseDates[i];
				regressionBasisfunctionTimes[numberOfPeriods] = swapEndDate;
				ArrayList<RandomVariable> basisFunctions = new ArrayList<>();

				double swapMaturity = swapEndDate;

				double exerciseTime = evaluationTime;

				int exerciseIndex = Arrays.binarySearch(regressionBasisfunctionTimes, exerciseTime);
				if (exerciseIndex < 0) {
					exerciseIndex = -exerciseIndex;
				}
				if (exerciseIndex >= numberOfPeriods) {
					exerciseIndex = numberOfPeriods - 1;
				}

				// Constant
				RandomVariable one = new RandomVariableFromDoubleArray(1.0);

				RandomVariable basisFunction = one;
				basisFunctions.add(basisFunction);

				// forward rate to the next period
				RandomVariable rateShort = model.getLIBOR(exerciseTime, exerciseTime,
						regressionBasisfunctionTimes[exerciseIndex + 1]);
				basisFunctions.add(rateShort);
				basisFunctions.add(rateShort.pow(2.0));

				// forward rate to the end of the product
				RandomVariable rateLong = model.getLIBOR(exerciseTime, regressionBasisfunctionTimes[exerciseIndex],
						swapMaturity);
				basisFunctions.add(rateLong);
				basisFunctions.add(rateLong.pow(2.0));

				// Numeraire (adapted to multicurve framework)
				RandomVariable discountFactor = model.getNumeraire(exerciseTime).invert();
				basisFunctions.add(discountFactor);

				// Cross
				basisFunctions.add(rateLong.mult(discountFactor));

				return basisFunctions.toArray(new RandomVariable[basisFunctions.size()]);
			}
		};
	}

	public RandomVariable getValueOfLegAnalytic(double evaluationTime, double[] exerciseDates, double[] paymentDates,
			LIBORModelMonteCarloSimulationModel model, boolean paysFloatingRate, double fixRate)
			throws CalculationException {

		RandomVariable discountedCashflowFloatingLeg = model.getRandomVariableForConstant(0.0);
		for (int periodIndex = exerciseDates.length - 1; periodIndex >= 0; periodIndex--) {
			double paymentTime = paymentDates[periodIndex];
			double fixingTime = exerciseDates[periodIndex];
			double periodLength = bermudanOption.getPeriodLengths()[periodIndex];
			double notional = bermudanOption.getPeriodNotionals()[periodIndex];
			/*
			 * Note that it is important that getForwardDiscountBond and getLIBOR are called
			 * with evaluationTime = exerciseTime.
			 */
			// RandomVariable discountBond =
			// model.getModel().getForwardDiscountBond(evaluationTime, paymentTime);
			if (paysFloatingRate) {
				RandomVariable libor = model.getLIBOR(evaluationTime, fixingTime, paymentTime);
				RandomVariable periodCashFlow = libor.mult(periodLength).mult(notional);
				discountedCashflowFloatingLeg = discountedCashflowFloatingLeg
						.add(periodCashFlow.div(model.getNumeraire(fixingTime)));
			}
			if (fixRate != 0) {
				RandomVariable periodCashFlow = model.getRandomVariableForConstant(fixRate * periodLength * notional);
				discountedCashflowFloatingLeg = discountedCashflowFloatingLeg
						.add(periodCashFlow.div(model.getNumeraire(fixingTime)));
			}
		}
		return discountedCashflowFloatingLeg;

	}

	// some getters
	@Override
	public RandomVariable[] getCacheValuesOfUnderlying() {
		return cacheValuesOfUnderlying;
	}

	@Override
	public RandomVariable[] getCacheConditionalExpectations() {
		return cacheConditionalExpectations;
	}

	public RegressionBasisFunctionsProvider getRegressionBasisFunctionsProvider() {
		return regressionBasisFunctionsProvider;
	}

	public void setRegressionBasisFunctionsProvider(RegressionBasisFunctionsProvider regressionBasisFunctionProvider) {
		this.regressionBasisFunctionsProvider = regressionBasisFunctionProvider;
	}
}