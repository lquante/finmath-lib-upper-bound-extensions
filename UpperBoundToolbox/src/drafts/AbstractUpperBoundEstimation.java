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

		// initialize cache arrays for calculation
		int numberOfPeriods = bermudanOption.getFixingDates().length;
		cacheOptionValues = new RandomVariable[numberOfPeriods];
		cacheValuesOfUnderlying = new RandomVariable[numberOfPeriods];
		cacheConditionalExpectations = new RandomVariable[numberOfPeriods];
		cacheTriggerValues = new RandomVariable[numberOfPeriods];
		optionValue = model.getRandomVariableForConstant(0);
		// calculate value of lower bound
		int evaluationTimeIndex = model.getTimeIndex(evaluationTime);

		this.bermudanOption.setValuationMethod(lowerBoundMethod);
		this.bermudanOption.getValue(evaluationTime, model);
		RandomVariable[] cacheUnderlying = ((SimpleLowerBoundEstimation) this.bermudanOption.getValuationMethod())
				.getCacheValuesOfUnderlying();
		RandomVariable[] cacheOptionValues = lowerBoundMethod.getCacheOptionValues();
		// Check exercise condition of lower bound method
		RandomVariable[] cacheTriggers = lowerBoundMethod.getCacheTriggerValues();
		// calculate upper bound

		double deltaZeroApproximation = calculateDeltaZero(evaluationTimeIndex, model, cacheUnderlying, cacheOptionValues,
				cacheTriggers);
		// Note that values is a relative price - no numeraire division is required
		RandomVariable numeraireAtEvaluationTime = model.getNumeraire(evaluationTime);
		RandomVariable monteCarloProbabilitiesAtEvaluationTime = model.getMonteCarloWeights(evaluationTime);
		RandomVariable discountFactor = numeraireAtEvaluationTime.div(monteCarloProbabilitiesAtEvaluationTime);

		optionValue = cacheOptionValues[evaluationTimeIndex].mult(discountFactor).add(deltaZeroApproximation);

		return optionValue;

	}

	protected abstract double calculateDeltaZero(int period, LIBORModelMonteCarloSimulationModel model,
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
