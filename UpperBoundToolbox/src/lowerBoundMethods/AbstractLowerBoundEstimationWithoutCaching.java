package lowerBoundMethods;

import bermudanSwaptionFramework.BermudanSwaption;
import bermudanSwaptionFramework.BermudanSwaptionValueEstimatorInterface;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.MonteCarloSimulationModel;
import net.finmath.montecarlo.automaticdifferentiation.RandomVariableDifferentiable;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAAD;
import net.finmath.montecarlo.conditionalexpectation.RegressionBasisFunctionsProvider;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import java.util.HashMap;
import java.util.Map;

/**
 * Provides a framework for a lower bound estimation without caching to increase performance, based on a specific method to calculate the exercise trigger values.
 * @author Lennart Quante
 * @version 1.0
 */
public abstract class AbstractLowerBoundEstimationWithoutCaching extends AbstractLowerBoundEstimation implements BermudanSwaptionValueEstimatorInterface {

	
	
	@Override
	public RandomVariable getValueEstimation(BermudanSwaption bermudanOption, double evaluationTime,
			LIBORModelMonteCarloSimulationModel model, RandomVariable triggerValuesInput) throws CalculationException {

		this.model = model;
		this.bermudanOption = bermudanOption;
		
		// initialize cache arrays for calculation
		int numberOfFixingDates = bermudanOption.getFixingDates().length - 1;
		
		initializeValuation (triggerValuesInput);
		return backwardAlgorithmValuation(evaluationTime);
		
	}

	@Override
	protected  RandomVariable calculateLiborRate(double fixingDate, double paymentDate) throws CalculationException {
		return  (model.getLIBOR(fixingDate, fixingDate,
				paymentDate));
	}
	@Override 
	void applyExerciseCriteria(double exerciseDate, int period)
	{
		optionValue = triggerValues.choose(exerciseValue, continuationValue).floor(0);

		exerciseTime = triggerValues.choose(new Scalar(exerciseDate), exerciseTime);
		
	}

	

	// abstract methods to be specified by the various implementations
	
	protected abstract RandomVariable calculateTriggerValues(int period, double fixingDate,
			LIBORModelMonteCarloSimulationModel model) throws CalculationException;

	// some getters

	

	public double[] getExerciseProbablities() {
		return exerciseProbablities;
	}

	public RandomVariable getExerciseTime() {
		return exerciseTime;
	}




}