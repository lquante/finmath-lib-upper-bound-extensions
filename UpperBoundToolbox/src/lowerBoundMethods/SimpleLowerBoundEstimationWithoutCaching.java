package lowerBoundMethods;

import java.util.ArrayList;
import java.util.Arrays;

import bermudanSwaptionFramework.SimplestExerciseStrategy;
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
 * Implements a lower bound estimation of a Bermudan Swaption value without caching.
 * @author (c)Christian P. Fries, modified by Lennart Quante
 *	
 */
public class SimpleLowerBoundEstimationWithoutCaching extends AbstractLowerBoundEstimationWithoutCaching
		implements RegressionBasisFunctionsProvider {

	public enum BasisFunctionType {
		SwapRates, ForwardRates
	}

	RandomVariable triggerValues;
	RegressionBasisFunctionsProvider regressionBasisFunctionsProvider;

	/**
	 * @param basisFunctionType choice of special basis functions to be used
	 */
	public SimpleLowerBoundEstimationWithoutCaching(BasisFunctionType basisFunctionType) {

		if (basisFunctionType == BasisFunctionType.ForwardRates)
			setRegressionBasisFunctionsProvider(getBasisFunctionsProviderWithSwapRates());
		if (basisFunctionType == BasisFunctionType.SwapRates)
			setRegressionBasisFunctionsProvider(getBasisFunctionsProviderWithForwardRates());

	}

	/**
	 * Generates value estimator with default basis functions
	 */
	public SimpleLowerBoundEstimationWithoutCaching() {

	}
	/* (non-Javadoc)
	 * @see drafts.AbstractSimpleBoundEstimation#calculateTriggerValues(int, double, net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel)
	 */
	@Override
	protected RandomVariable calculateTriggerValues(int period, double fixingDate,
			LIBORModelMonteCarloSimulationModel model) throws CalculationException {

		SimplestExerciseStrategy exerciseStrategy = new SimplestExerciseStrategy();
		// Calculate the exercise criteria (exercise if the following trigger is
		// negative)
		RandomVariable triggerValuesDiscounted = exerciseStrategy.getTriggerValues (continuationValue,
				exerciseValue);

		// Remove foresight through condition expectation
		ConditionalExpectationEstimator conditionalExpectationOperator = getConditionalExpectationEstimator(fixingDate,
				model);

		// Calculate conditional expectation. Note that no discounting (numeraire division) is required!
		triggerValues = triggerValuesDiscounted.getConditionalExpectation(conditionalExpectationOperator);
		
		return triggerValues;

	}

	

	// some getters and setters
	
	public RandomVariable getTriggerValues() {
		return triggerValues;
	}

	public void setTriggerValues(RandomVariableDifferentiableAAD triggerValues) {
		this.triggerValues = triggerValues;
	}
	

}