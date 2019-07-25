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
	RandomVariable[] cacheConditionalExpectations;

	RandomVariable continuationValue;
	RandomVariable exerciseValue;
	RandomVariable optionValue;
	RandomVariable triggerValues;

	private AbstractSimpleBoundEstimation lowerBoundMethod;
	RandomVariable[] cacheTriggerValues;

	public AbstractUpperBoundEstimation(AbstractSimpleBoundEstimation lowerBoundMethod) {

		this.lowerBoundMethod = lowerBoundMethod;
	}
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
			LIBORModelMonteCarloSimulationModel model, RandomVariable triggerValuesInput) throws CalculationException {

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
		cacheTriggerValues = new RandomVariable[numberOfPeriods + 1];
		optionValue = model.getRandomVariableForConstant(0);
		// calculate value of lower bound
		int period = 0;
		double valuationTime = this.bermudanOption.getFixingDates()[period];
		this.bermudanOption.setValuationMethod(lowerBoundMethod);
		this.bermudanOption.getValue(valuationTime, model);
		RandomVariable[] cacheUnderlying = ((SimpleLowerBoundEstimation) this.bermudanOption.getValuationMethod())
				.getCacheValuesOfUnderlying();
		RandomVariable[] cacheOptionValues = lowerBoundMethod.getCacheOptionValues();
		// Check exercise condition of lower bound method
		RandomVariable[] cacheTriggers = lowerBoundMethod.getCacheTriggerValues();
		// calculate upper bound
		if (this.bermudanOption.getIsPeriodStartDateExerciseDate()[period])
			optionValue = calculateOptionValue(period, model, cacheUnderlying, cacheOptionValues, cacheTriggers);
		// the returned upper bound has to be already discounted to evaluation time and
		// monte carlo weights applied
		return optionValue;
	}

	protected abstract RandomVariable calculateOptionValue(int period, LIBORModelMonteCarloSimulationModel model,
			RandomVariable[] cacheUnderlying, RandomVariable[] cacheOptionValues, RandomVariable[] triggerValues)
			throws CalculationException;

	public RandomVariable[] getCacheValuesOfUnderlying() {
		return (cacheValuesOfUnderlying);
	}

	public RandomVariable[] getCacheConditionalExpectations() {
		return (cacheConditionalExpectations);
	}

	public RandomVariable[] getCacheOptionValues() {
		return (cacheOptionValues);
	}

	public RandomVariable getTriggerValues() {
		return triggerValues;
	}

}
