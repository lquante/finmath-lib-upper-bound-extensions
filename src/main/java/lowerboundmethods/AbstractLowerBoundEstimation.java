package lowerboundmethods;

import java.util.ArrayList;
import java.util.Arrays;

import bermudanswaptionframework.BermudanSwaption;
import bermudanswaptionframework.BermudanSwaptionValueEstimatorInterface;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.MonteCarloSimulationModel;
import net.finmath.montecarlo.RandomVariableFromDoubleArray;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAAD;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import net.finmath.montecarlo.conditionalexpectation.RegressionBasisFunctionsProvider;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.ConditionalExpectationEstimator;
import net.finmath.stochastic.RandomVariable;

/**
 * Provides a framework for a lower bound estimation via the backward algorithm,
 * based on a specific method to calculate the exercise trigger values.
 * 
 * @author Lennart Quante
 * @version 1.0
 */
public abstract class AbstractLowerBoundEstimation implements BermudanSwaptionValueEstimatorInterface {

	RegressionBasisFunctionsProvider regressionBasisFunctionsProvider;
	BermudanSwaption bermudanSwaption;

	RandomVariable continuationValue;
	RandomVariable exerciseValue;
	RandomVariable optionValue;
	RandomVariable triggerValues;


	RandomVariable exerciseTime;
	protected LIBORModelMonteCarloSimulationModel model;

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * drafts.BermudanSwaptionValueEstimatorInterface#getValueEstimation(drafts.
	 * BermudanSwaption, double,
	 * net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel)
	 */

	// general methods to avoid code duplication in classes inheriting this.
	protected void initializeValuation(RandomVariable triggerValuesInput) {

		// After the last period the product has value zero: Initialize values to zero.
		RandomVariable zero = model.getRandomVariableForConstant(0.0);
		exerciseValue = zero;
		continuationValue = zero;

		exerciseTime = model.getRandomVariableForConstant(Double.POSITIVE_INFINITY);
		// use input trigger values if given
		triggerValues = triggerValuesInput;

	}

	/**
	 * Valuation of the Bermudan Swaption via the backward valuation algorithm.
	 * @param evaluationTime The time at which the backwardAlgorthm should evaluate the option.
	 * @return The value random variable
	 * @throws CalculationException
	 */
	protected RandomVariable backwardAlgorithmValuation(double evaluationTime) throws CalculationException {
		int numberOfFixingDates = bermudanSwaption.getFixingDates().length - 1;
		// Loop backward over the swap periods
		for (int period = numberOfFixingDates; period >= 0; period--) {
			double fixingDate = bermudanSwaption.getFixingDates()[period];
			double exerciseDate = fixingDate;
			double periodLength = bermudanSwaption.getPeriodLengths()[period];
			double paymentDate = bermudanSwaption.getPaymentDates()[period];
			double notional = bermudanSwaption.getPeriodNotionals()[period];
			double swaprate = bermudanSwaption.getSwaprates()[period];
			// Get random variables - note that this is the rate at simulation time = exerciseDate
			RandomVariable libor = calculateLiborRate(fixingDate, paymentDate);
			// calculate payoff
			RandomVariable payoff = libor.sub(swaprate).mult(periodLength).mult(notional);
			// Apply discounting and Monte-Carlo probabilities
			RandomVariable numeraire = model.getNumeraire(paymentDate);
			RandomVariable monteCarloProbabilities = model.getMonteCarloWeights(paymentDate);
			payoff = payoff.div(numeraire).mult(monteCarloProbabilities);
			if (bermudanSwaption.isCallable()) 
				exerciseValue = exerciseValue.add(payoff);
			else
				continuationValue = continuationValue.add(payoff); // cancelable
			determineTriggerValuesWithEarlyTermination(period,fixingDate);
			applyExerciseCriteria(exerciseDate, period);
		}
		// Note that value has to be a relative price - no numeraire division is performed!
		RandomVariable numeraireAtZero = model.getNumeraire(evaluationTime);
		RandomVariable monteCarloProbabilitiesAtZero = model.getMonteCarloWeights(evaluationTime);
		return optionValue.mult(numeraireAtZero).div(monteCarloProbabilitiesAtZero);
	}
	// method to maybe enhance techniques of early termination
	private void determineTriggerValuesWithEarlyTermination(int period, double fixingDate) throws CalculationException {
		if (this.bermudanSwaption.getIsPeriodStartDateExerciseDate()[period]) {
			if(bermudanSwaption.isCallable()==false)
				triggerValues = calculateTriggerValues(period, fixingDate, model);
			else
				// small extension to account for points where exercise won't occur (compare Joshi2014)
				if (exerciseValue.getMax()<=0)
					triggerValues= model.getRandomVariableForConstant(-1);
				else
					triggerValues = calculateTriggerValues(period, fixingDate, model);
		}
	}
	// abstract methods to be specified by the various implementations
	@Override
	public abstract RandomVariable getValueEstimation(BermudanSwaption bermudanOption, double evaluationTime,
			LIBORModelMonteCarloSimulationModel model, RandomVariable triggerValuesInput) throws CalculationException;

	/**
	 * Method to retrieve the LIBOR rate for the input parameters.
	 * @param fixingDate
	 * @param paymentDate
	 * @return The LIBOR rate from fixing date to payment date
	 * @throws CalculationException
	 */
	protected abstract RandomVariable calculateLiborRate(double fixingDate, double paymentDate)
			throws CalculationException;

	/**
	 * Returns the triggerValues for the input paramters using the exerciseStrategy of the valuation interface, to be specified by specific method.
	 * @param period
	 * @param fixingDate
	 * @param model
	 * @return RandomVariable of triggerValues
	 * @throws CalculationException
	 */
	protected abstract RandomVariable calculateTriggerValues(int period, double fixingDate,
			LIBORModelMonteCarloSimulationModel model) throws CalculationException;

	/**
	 * Applies the exercise criteria to the option.
	 * @param exerciseDate
	 * @param period
	 */
	abstract void applyExerciseCriteria(double exerciseDate, int period);

	// basis functions used in the algorithms
	/**
	 * Return the conditional expectation estimator suitable for this product.
	 *
	 * @param fixingDate The condition time.
	 * @param model      The model
	 * @return The conditional expectation estimator suitable for this product
	 * @throws net.finmath.exception.CalculationException Thrown if the valuation
	 *                                                    fails, specific cause may
	 *                                                    be available via the
	 *                                                    <code>cause()</code>
	 *                                                    method.
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see drafts.AbstractSimpleBoundEstimation#getBasisFunctions(double,
	 * net.finmath.montecarlo.MonteCarloSimulationModel)
	 */

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
	public RandomVariable[] getBasisFunctions(double evaluationTime, LIBORModelMonteCarloSimulationModel model)
			throws CalculationException {
		double[] fixingDates = bermudanSwaption.getFixingDates();
		double[] paymentDates = bermudanSwaption.getPaymentDates();
		ArrayList<RandomVariable> basisFunctions = standardBasisFunctions(evaluationTime, model,fixingDates,paymentDates);
		return basisFunctions.toArray(new RandomVariable[basisFunctions.size()]);
	}


	/**
	 * Provides basis functions based on a constant 1, a rate to the next period, a long rate and the numeraire
	 * @param evaluationTime
	 * @param model
	 * @param fixingDates
	 * @param paymentDates
	 * @return ArrayList<RandomVariable> of the estimated basis functions.
	 * @throws CalculationException
	 */
	private ArrayList<RandomVariable> standardBasisFunctions(double evaluationTime, LIBORModelMonteCarloSimulationModel model,double[] fixingDates, double[] paymentDates)
			throws CalculationException {
		ArrayList<RandomVariable> basisFunctions = new ArrayList<>();
		// Constant
		RandomVariable one = new RandomVariableFromDoubleArray(1.0);
		basisFunctions.add(one);
		int fixingDateIndex = Arrays.binarySearch(fixingDates, evaluationTime);
		if (fixingDateIndex < 0)
			fixingDateIndex = -fixingDateIndex;
		if (fixingDateIndex >= fixingDates.length)
			fixingDateIndex = fixingDates.length - 1;
		// forward rate to the next period
		RandomVariable rateShort = model.getLIBOR(evaluationTime, evaluationTime, paymentDates[fixingDateIndex]);
		RandomVariable discountShort = rateShort.mult(paymentDates[fixingDateIndex] - evaluationTime).add(1.0).invert();
		basisFunctions.add(discountShort);
		basisFunctions.add(discountShort.pow(2.0));
		// forward rate to the end of the product
		RandomVariable rateLong = model.getLIBOR(evaluationTime, fixingDates[fixingDateIndex],
				paymentDates[paymentDates.length - 1]);
		RandomVariable discountLong = rateLong
				.mult(paymentDates[paymentDates.length - 1] - fixingDates[fixingDateIndex]).add(1.0).invert();
		basisFunctions.add(discountLong);
		basisFunctions.add(discountLong.pow(2.0));
		// Numeraire
		RandomVariable numeraire = model.getNumeraire(evaluationTime).invert();
		basisFunctions.add(numeraire);
		return basisFunctions;
	}

	/*
	 * Some popular variants to create regression basis functions
	 */

	/**
	 *
	 * @return An \( \mathcal{F}_{t} \)-measurable random variable, basis functions
	 *         using swap rate estimation
	 * 
	 */
	public RegressionBasisFunctionsProvider getBasisFunctionsProviderWithSwapRates() {
		return new RegressionBasisFunctionsProvider() {
			@Override
			public RandomVariableDifferentiableAAD[] getBasisFunctions(double evaluationTime,
					MonteCarloSimulationModel monteCarloModel) throws CalculationException {
				LIBORModelMonteCarloSimulationModel model = (LIBORModelMonteCarloSimulationModel) monteCarloModel;

				double[] fixingDates = bermudanSwaption.getFixingDates();
				double[] regressionBasisfunctionTimes = bermudanSwaption.getPaymentDates();

				ArrayList<RandomVariable> basisFunctions = new ArrayList<>();
				basisFunctions.addAll(standardBasisFunctions(evaluationTime, model, fixingDates, regressionBasisfunctionTimes));
				
				basisFunctions.addAll(swapRateBasisFunctions(evaluationTime, model, fixingDates, regressionBasisfunctionTimes));
				
				return basisFunctions.toArray(new RandomVariableDifferentiableAAD[basisFunctions.size()]);
			}
		};
	}
	
	/**
	 * Provides basis functions based on a swap rate estimate.
	 * @param evaluationTime
	 * @param model
	 * @param fixingDates
	 * @param regressionBasisfunctionTimes
	 * @return ArrayList<RandomVariable> of the estimated basis functions.
	 * @throws CalculationException
	 */
	private ArrayList<RandomVariable> swapRateBasisFunctions(double evaluationTime, LIBORModelMonteCarloSimulationModel model,double[] fixingDates, double[] regressionBasisfunctionTimes)
			throws CalculationException {
		ArrayList<RandomVariable> basisFunctions = new ArrayList<>();
		double exerciseTime = evaluationTime;
		int exerciseIndex = Arrays.binarySearch(regressionBasisfunctionTimes, exerciseTime);
		if (exerciseIndex < 0) 
			exerciseIndex = -exerciseIndex;
		if (exerciseIndex >= fixingDates.length) 
			exerciseIndex = fixingDates.length - 1;
		// Numeraire (adapted to multicurve framework)
		RandomVariable discountFactor = model.getNumeraire(exerciseTime).invert();
		// Add swap rates of underlyings.
		for (int exerciseIndexUnderlying = exerciseIndex; exerciseIndexUnderlying < fixingDates.length; exerciseIndexUnderlying++) {
			RandomVariable floatLeg = getValueOfLegAnalytic(exerciseTime, fixingDates, regressionBasisfunctionTimes, model, true, 0.0);
			RandomVariable annuity = getValueOfLegAnalytic(exerciseTime, fixingDates, regressionBasisfunctionTimes, model, false, 1.0);
			RandomVariable swapRate = floatLeg.div(annuity);
			RandomVariable basisFunction = swapRate.mult(discountFactor);
			basisFunctions.add(basisFunction);
			basisFunctions.add(basisFunction.squared());
		}
		return basisFunctions;
	}

	/**
	 * 
	 * @return An \( \mathcal{F}_{t} \)-measurable random variable, basis functions
	 *         using forward rates and a cross.
	 * 
	 */
	public RegressionBasisFunctionsProvider getBasisFunctionsProviderWithForwardRates() {
		return new RegressionBasisFunctionsProvider() {
			@Override
			public RandomVariable[] getBasisFunctions(double evaluationTime, MonteCarloSimulationModel monteCarloModel) throws CalculationException {
				LIBORModelMonteCarloSimulationModel model = (LIBORModelMonteCarloSimulationModel) monteCarloModel;
				double[] fixingDates = bermudanSwaption.getFixingDates();
				double swapEndDate = bermudanSwaption.getPaymentDates()[fixingDates.length - 1];
				double[] regressionBasisfunctionTimes = new double[fixingDates.length + 1];
				int numberOfPeriods = fixingDates.length;
				for (int i = 0; i < fixingDates.length; i++)
					regressionBasisfunctionTimes[i] = fixingDates[i];
				regressionBasisfunctionTimes[numberOfPeriods] = swapEndDate;
				ArrayList<RandomVariable> basisFunctions = new ArrayList<>();
				basisFunctions.addAll(standardBasisFunctions(evaluationTime, model, fixingDates, regressionBasisfunctionTimes));
				basisFunctions.addAll(crossBasisFunctions(swapEndDate,evaluationTime,model,regressionBasisfunctionTimes,numberOfPeriods));
				return basisFunctions.toArray(new RandomVariableDifferentiableAAD[basisFunctions.size()]);			
			}	
		};
	}

	/**
	 * Provides basis functions with a cross of the numeraire and a long rate.
	 * @param evaluationTime
	 * @param model
	 * @param fixingDates
	 * @param regressionBasisfunctionTimes
	 * @return ArrayList of RandomVariable of the estimated basis functions.
	 * @throws CalculationException
	 */
	private ArrayList<RandomVariable> crossBasisFunctions(double swapEndDate, double evaluationTime,
			LIBORModelMonteCarloSimulationModel model, double[] regressionBasisfunctionTimes, int numberOfPeriods) throws CalculationException {
		double swapMaturity = swapEndDate;
		double exerciseTime = evaluationTime;
		ArrayList<RandomVariable> basisFunctions = new ArrayList<>();
		int exerciseIndex = Arrays.binarySearch(regressionBasisfunctionTimes, exerciseTime);
		if (exerciseIndex < 0) 
			exerciseIndex = -exerciseIndex;
		if (exerciseIndex >= numberOfPeriods) 
			exerciseIndex = numberOfPeriods - 1;
		// forward rate to the end of the product
		RandomVariable rateLong = model.getLIBOR(exerciseTime, regressionBasisfunctionTimes[exerciseIndex],swapMaturity);
		basisFunctions.add(rateLong);
		basisFunctions.add(rateLong.pow(2.0));
		// Numeraire (adapted to multicurve framework)
		RandomVariable discountFactor = model.getNumeraire(exerciseTime).invert();
		basisFunctions.add(discountFactor);
		basisFunctions.add(rateLong.mult(discountFactor)); // Cross
		return basisFunctions;
	}
	
	/**
	 * Helper method to calculate analytic leg value for swap rate estimate.
	 * @param evaluationTime
	 * @param exerciseDates
	 * @param paymentDates
	 * @param model
	 * @param paysFloatingRate
	 * @param fixRate
	 * @return The value of the floating leg.
	 * @throws CalculationException
	 */
	public RandomVariable getValueOfLegAnalytic(double evaluationTime, double[] exerciseDates, double[] paymentDates,
			LIBORModelMonteCarloSimulationModel model, boolean paysFloatingRate, double fixRate)
					throws CalculationException {

		RandomVariable discountedCashflowFloatingLeg = model.getRandomVariableForConstant(0.0);
		for (int periodIndex = exerciseDates.length - 1; periodIndex >= 0; periodIndex--) {
			double paymentTime = paymentDates[periodIndex];
			double fixingTime = exerciseDates[periodIndex];
			double periodLength = bermudanSwaption.getPeriodLengths()[periodIndex];
			double notional = bermudanSwaption.getPeriodNotionals()[periodIndex];
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


	public RandomVariable getExerciseTime() {
		return exerciseTime;
	}

	public RegressionBasisFunctionsProvider getRegressionBasisFunctionsProvider() {
		return regressionBasisFunctionsProvider;
	}

	public void setRegressionBasisFunctionsProvider(RegressionBasisFunctionsProvider regressionBasisFunctionProvider) {
		this.regressionBasisFunctionsProvider = regressionBasisFunctionProvider;
	}

}