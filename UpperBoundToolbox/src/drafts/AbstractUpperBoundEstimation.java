package drafts;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.conditionalexpectation.RegressionBasisFunctionsProvider;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;

public abstract class AbstractUpperBoundEstimation implements BermudanSwaptionValueEstimatorInterface {
	// basis functions used in the lower bound
	protected RegressionBasisFunctionsProvider regressionBasisFunctionsProvider;
	// product to be evaluated
	BermudanSwaption bermudanOption;
	// cache arrays for upper bound estimation
	RandomVariable[] cacheOptionValues;
	RandomVariable[] cacheValuesOfUnderlying;

	public RandomVariable getTriggerValues() {
		return triggerValues;
	}

	RandomVariable continuationValue;
	RandomVariable exerciseValue;
	RandomVariable optionValue;
	RandomVariable triggerValues;

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * drafts.BermudanSwaptionValueEstimatorInterface#getValueEstimation(drafts.
	 * BermudanSwaption, double,
	 * net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel)
	 */
	@Override
	public RandomVariable getValueEstimation(BermudanSwaption bermudanOption, double evaluationTime,
			LIBORModelMonteCarloSimulationModel model) throws CalculationException {

		this.bermudanOption = bermudanOption;

		// After the last period the product has value zero: Initialize values to zero.
		exerciseValue = model.getRandomVariableForConstant(0.0);
		continuationValue = model.getRandomVariableForConstant(0.0);

		RandomVariable exerciseTime = model.getRandomVariableForConstant(Double.POSITIVE_INFINITY);
		triggerValues = null;

		// initialize cache arrays for calculation
		int numberOfPeriods = bermudanOption.getFixingDates().length - 1;
		cacheOptionValues = new RandomVariable[numberOfPeriods + 1];
		cacheValuesOfUnderlying = new RandomVariable[numberOfPeriods + 1];
		cacheConditionalExpectations = new RandomVariable[numberOfPeriods + 1];
		// Loop backward over the swap periods
		for (int period = numberOfPeriods; period >= 0; period--) {

			if (this.bermudanOption.getIsPeriodStartDateExerciseDate()[period]) {
				optionValue = calculateOptionValue(period, model);
			}
			cacheOptionValues[period] = optionValue;

		}
		// the returned upper bound has to be already discounted to evaluation time and
		// monte carlo weights are applied
		return optionValue;
	}

	protected abstract RandomVariable calculateOptionValue(int period, LIBORModelMonteCarloSimulationModel model)
			throws CalculationException;

	public RandomVariable[] getCacheValuesOfUnderlying() {
		return cacheValuesOfUnderlying;
	}

	public RandomVariable[] getCacheConditionalExpectations() {
		return cacheConditionalExpectations;
	}

	RandomVariable[] cacheConditionalExpectations;

	public RandomVariable[] getCacheOptionValues() {
		return cacheOptionValues;
	}
}
