package drafts;

import java.util.HashSet;
import java.util.Set;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.MonteCarloSimulationModel;
import net.finmath.montecarlo.automaticdifferentiation.RandomVariableDifferentiable;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAAD;
import net.finmath.montecarlo.conditionalexpectation.RegressionBasisFunctionsProvider;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

public abstract class AbstractSimpleBoundEstimation implements BermudanSwaptionValueEstimatorInterface {

	RegressionBasisFunctionsProvider regressionBasisFunctionsProvider;
	BermudanSwaption bermudanOption;
	RandomVariable[] cacheOptionValues;
	RandomVariable[] cacheValuesOfUnderlying;
	RandomVariable[] cacheConditionalExpectations;
	Set<Long> liborIDs;

	RandomVariable continuationValue;
	RandomVariable exerciseValue;
	RandomVariable optionValue;
	RandomVariable triggerValues;
	RandomVariable[] cacheTriggerValues;
	double[] exerciseProbablities;
	RandomVariable exerciseTime;

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
		liborIDs = new HashSet<>();
		// After the last period the product has value zero: Initialize values to zero.
		RandomVariable zero = model.getRandomVariableForConstant(0.0);
		exerciseValue = zero;
		continuationValue = zero;

		exerciseTime = model.getRandomVariableForConstant(Double.POSITIVE_INFINITY);
		triggerValues = triggerValuesInput;

		// initialize cache arrays for calculation
		int numberOfPeriods = bermudanOption.getFixingDates().length - 1;
		cacheOptionValues = new RandomVariable[numberOfPeriods + 1];
		cacheValuesOfUnderlying = new RandomVariable[numberOfPeriods + 1];
		cacheConditionalExpectations = new RandomVariable[numberOfPeriods + 1];
		cacheTriggerValues = new RandomVariable[numberOfPeriods + 1];

		exerciseProbablities = new double[numberOfPeriods + 1];
		// Loop backward over the swap periods
		for (int period = numberOfPeriods; period >= 0; period--) {
			double fixingDate = bermudanOption.getFixingDates()[period];
			double exerciseDate = fixingDate;
			double periodLength = bermudanOption.getPeriodLengths()[period];
			double paymentDate = bermudanOption.getPaymentDates()[period];
			double notional = bermudanOption.getPeriodNotionals()[period];
			double swaprate = bermudanOption.getSwaprates()[period];

			// Get random variables - note that this is the rate at simulation time =
			// exerciseDate
			RandomVariable libor = ((RandomVariableDifferentiableAAD) model.getLIBOR(fixingDate, fixingDate,
					fixingDate + periodLength)).getCloneIndependent();
			Long liborID = ((RandomVariableDifferentiable) libor).getID();
			liborIDs.add(liborID);
			RandomVariable payoff = libor.sub(swaprate).mult(periodLength).mult(notional);

			// Apply discounting and Monte-Carlo probabilities
			RandomVariable numeraire = model.getNumeraire(paymentDate);
			RandomVariable monteCarloProbabilities = model.getMonteCarloWeights(paymentDate);
			payoff = payoff.div(numeraire).mult(monteCarloProbabilities);

			if (bermudanOption.isCallable()) {
				exerciseValue = exerciseValue.add(payoff);
			} else {
				continuationValue = continuationValue.add(payoff); // cancelable
			}
			// Calculate the exercise criteria (exercise if the following trigger is
			// negative)
			if (this.bermudanOption.getIsPeriodStartDateExerciseDate()[period]) {
				triggerValues = calculateTriggerValues(period, fixingDate, model);
			}
			// Apply the exercise criteria

			optionValue = triggerValues.choose(exerciseValue, continuationValue).floor(0);

			exerciseTime = triggerValues.choose(exerciseTime, new Scalar(exerciseDate));

			// caching for upper bound methods
			cacheTriggerValues[period] = triggerValues;
			cacheOptionValues[period] = optionValue;

		}

		// exerciseProbablities =
		// exerciseTime.getHistogram(bermudanOption.getFixingDates());

		// Note that values is a relative price - no numeraire division is required
		RandomVariable numeraireAtZero = model.getNumeraire(evaluationTime);
		RandomVariable monteCarloProbabilitiesAtZero = model.getMonteCarloWeights(evaluationTime);

		optionValue = optionValue.mult(numeraireAtZero).div(monteCarloProbabilitiesAtZero);

		return optionValue;
	}

	protected abstract RandomVariable calculateTriggerValues(int period, double fixingDate,
			LIBORModelMonteCarloSimulationModel model) throws CalculationException;

	// some getters

	public abstract RandomVariable[] getBasisFunctions(double evaluationTime, MonteCarloSimulationModel model)
			throws CalculationException;

	public double[] getExerciseProbablities() {
		return exerciseProbablities;
	}

	public RandomVariable getExerciseTime() {
		return exerciseTime;
	}

	public Set<Long> getLiborIDs() {
		return liborIDs;
	}

	public RandomVariable[] getCacheValuesOfUnderlying() {
		return cacheValuesOfUnderlying;
	}

	public RandomVariable[] getCacheConditionalExpectations() {
		return cacheConditionalExpectations;
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