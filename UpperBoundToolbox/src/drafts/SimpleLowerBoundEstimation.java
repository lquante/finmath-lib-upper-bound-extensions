package drafts;

import java.util.ArrayList;
import java.util.Arrays;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariableFromDoubleArray;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.ConditionalExpectationEstimator;
import net.finmath.stochastic.RandomVariable;

public class SimpleLowerBoundEstimation extends AbstractSimpleBoundEstimation {

	RandomVariable triggerValues;

	@Override
	public RandomVariable getTriggerValues() {
		return triggerValues;
	}

	public void setTriggerValues(RandomVariable triggerValues) {
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

	/**
	 * Return the basis functions for the regression suitable for this product.
	 *
	 * @param fixingDate     The condition time.
	 * @param model          The model
	 * @param bermudanOption
	 * @return The basis functions for the regression suitable for this product.
	 * @throws net.finmath.exception.CalculationException Thrown if the valuation
	 *         fails, specific cause may be available via the <code>cause()</code>
	 *         method.
	 */
	public RandomVariable[] getBasisFunctions(double fixingDate, LIBORModelMonteCarloSimulationModel model)
			throws CalculationException {

		ArrayList<RandomVariable> basisFunctions = new ArrayList<>();
		// Constant
		RandomVariable basisFunction = new RandomVariableFromDoubleArray(1.0);
		basisFunctions.add(basisFunction);

		int fixingDateIndex = Arrays.binarySearch(bermudanOption.getFixingDates(), fixingDate);
		if (fixingDateIndex < 0) {
			fixingDateIndex = -fixingDateIndex;
		}
		if (fixingDateIndex >= bermudanOption.getFixingDates().length) {
			fixingDateIndex = bermudanOption.getFixingDates().length - 1;
		}

		// forward rate to the next period
		RandomVariable rateShort = model.getLIBOR(fixingDate, fixingDate,
				bermudanOption.getPaymentDates()[fixingDateIndex]);
		RandomVariable discountShort = rateShort.mult(bermudanOption.getPaymentDates()[fixingDateIndex] - fixingDate)
				.add(1.0).invert();
		basisFunctions.add(discountShort);
		basisFunctions.add(discountShort.pow(2.0));
		// basisFunctions.add(rateShort.pow(3.0));

		// forward rate to the end of the product
		RandomVariable rateLong = model.getLIBOR(fixingDate, bermudanOption.getFixingDates()[fixingDateIndex],
				bermudanOption.getPaymentDates()[bermudanOption.getPaymentDates().length - 1]);
		RandomVariable discountLong = rateLong
				.mult(bermudanOption.getPaymentDates()[bermudanOption.getPaymentDates().length - 1]
						- bermudanOption.getFixingDates()[fixingDateIndex])
				.add(1.0).invert();
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

	@Override
	public RandomVariable[] getCacheValuesOfUnderlying() {
		return cacheValuesOfUnderlying;
	}

	@Override
	public RandomVariable[] getCacheConditionalExpectations() {
		return cacheConditionalExpectations;
	}
}