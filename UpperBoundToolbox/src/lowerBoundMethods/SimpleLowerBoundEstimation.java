package lowerBoundMethods;

import bermudanSwaptionFramework.SimplestExerciseStrategy;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAAD;
import net.finmath.montecarlo.conditionalexpectation.RegressionBasisFunctionsProvider;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.ConditionalExpectationEstimator;
import net.finmath.stochastic.RandomVariable;

/**
 * Implements a lower bound estimation of a Bermudan Swaption value which caches
 * values for reuse in upper bound methods.
 * 
 * @author (c)Christian P. Fries, modified by Lennart Quante
 * 
 */
public class SimpleLowerBoundEstimation extends AbstractLowerBoundEstimationInputForUpperBound
		implements RegressionBasisFunctionsProvider {

	public enum BasisFunctionType {
		SwapRates, ForwardRates
	}

	RandomVariable triggerValues;
	RegressionBasisFunctionsProvider regressionBasisFunctionsProvider;

	/**
	 * @param basisFunctionType choice of special basis functions to be used
	 */
	public SimpleLowerBoundEstimation(BasisFunctionType basisFunctionType) {

		if (basisFunctionType == BasisFunctionType.ForwardRates)
			setRegressionBasisFunctionsProvider(getBasisFunctionsProviderWithSwapRates());
		if (basisFunctionType == BasisFunctionType.SwapRates)
			setRegressionBasisFunctionsProvider(getBasisFunctionsProviderWithForwardRates());

	}

	/**
	 * Generates value estimator with default basis functions
	 */
	public SimpleLowerBoundEstimation() {

	}

	/**
	 * @see lowerBoundMethods.AbstractLowerBoundEstimation#calculateTriggerValues(int, double,
	 * net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel)
	 */
	@Override
	protected RandomVariable calculateTriggerValues(int period, double fixingDate,
			LIBORModelMonteCarloSimulationModel model) throws CalculationException {

		SimplestExerciseStrategy exerciseStrategy = new SimplestExerciseStrategy();
		// Calculate the exercise criteria (exercise if the following trigger is
		// negative)
		RandomVariable triggerValuesDiscounted = exerciseStrategy.getTriggerValues(continuationValue, exerciseValue);

		// Remove foresight through condition expectation
		ConditionalExpectationEstimator conditionalExpectationOperator = getConditionalExpectationEstimator(fixingDate,
				model);

		// Calculate conditional expectation. Note that no discounting (numeraire
		// division) is required!
		triggerValues = triggerValuesDiscounted.getConditionalExpectation(conditionalExpectationOperator);

		// cache values of underlying (including perfect foresight) and E[U_i|F_T_i-1]

		cacheValuesOfUnderlying[period] = exerciseValue;
		return triggerValues;

	}

	// some getters and setters

	@Override
	public RandomVariable getTriggerValues() {
		return triggerValues;
	}

	public void setTriggerValues(RandomVariableDifferentiableAAD triggerValues) {
		this.triggerValues = triggerValues;
	}

	@Override
	public RandomVariable[] getCacheValuesOfUnderlying() {
		return cacheValuesOfUnderlying;
	}

}