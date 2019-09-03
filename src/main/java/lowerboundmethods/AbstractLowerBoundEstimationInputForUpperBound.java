package lowerboundmethods;

import java.util.HashMap;
import java.util.Map;

import bermudanswaptionframework.BermudanSwaption;
import bermudanswaptionframework.BermudanSwaptionValueEstimatorInterface;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.automaticdifferentiation.RandomVariableDifferentiable;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAAD;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

/**
 * Provides a framework for a lower bound estimation with caching used as input
 * for upper bound methods, based on a specific method to calculate the exercise
 * trigger values.
 * 
 * @author Lennart Quante
 * @version 1.0
 */
public abstract class AbstractLowerBoundEstimationInputForUpperBound extends AbstractLowerBoundEstimation
		implements BermudanSwaptionValueEstimatorInterface {

	RandomVariable[] cacheOptionValues;
	RandomVariable[] cacheValuesOfUnderlying;
	Map<Double, Long> liborIDs;
	RandomVariable[] cacheTriggerValues;

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
		this.model = model;
		this.bermudanSwaption = bermudanOption;
		liborIDs = new HashMap<>();
		// initialize cache arrays for calculation
		int numberOfFixingDates = bermudanOption.getFixingDates().length - 1;
		cacheOptionValues = new RandomVariable[numberOfFixingDates + 1];
		cacheValuesOfUnderlying = new RandomVariable[numberOfFixingDates + 1];
		cacheTriggerValues = new RandomVariable[numberOfFixingDates + 1];
		initializeValuation(triggerValuesInput);
		optionValue= backwardAlgorithmValuation(evaluationTime);
	//	this.bermudanSwaption.setExerciseProbabilities(exerciseTime.getHistogram(this.bermudanSwaption.getFixingDates()));
		return optionValue;

	}

	@Override
	protected RandomVariable calculateLiborRate(double fixingDate, double paymentDate) throws CalculationException {
		RandomVariable libor = ((RandomVariableDifferentiableAAD) model.getLIBOR(fixingDate, fixingDate, paymentDate))
				.getCloneIndependent();
		// store liborIDs for use in automatic differentiation methods
		Long liborID = ((RandomVariableDifferentiable) libor).getID();
		liborIDs.put(fixingDate, liborID);
		return libor;
	}

	@Override
	void applyExerciseCriteria(double exerciseDate, int period) {
		optionValue = triggerValues.choose(exerciseValue, continuationValue).floor(0);

		exerciseTime = triggerValues.choose(new Scalar(exerciseDate), exerciseTime);

		// caching for upper bound methods
		cacheValuesOfUnderlying[period] = exerciseValue;
		cacheTriggerValues[period] = triggerValues;
		cacheOptionValues[period] = optionValue;
	}

	public Map<Double, Long> getLiborIDs() {
		return liborIDs;
	}

	public RandomVariable[] getCacheValuesOfUnderlying() {
		return cacheValuesOfUnderlying;
	}

	public RandomVariable[] getCacheOptionValues() {
		return cacheOptionValues;
	}

	public RandomVariable getTriggerValues() {
		return triggerValues;
	}

	public RandomVariable[] getCacheTriggerValues() {
		return cacheTriggerValues;
	}
}