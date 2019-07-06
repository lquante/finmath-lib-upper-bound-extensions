package drafts;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.MonteCarloSimulationModel;
import net.finmath.montecarlo.conditionalexpectation.RegressionBasisFunctionsProvider;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

public abstract class AbstractSimpleBoundEstimation implements BermudanSwaptionValueEstimatorInterface {

	RegressionBasisFunctionsProvider regressionBasisFunctionsProvider;
	BermudanSwaption bermudanOption;
	RandomVariable[] cacheOptionValues;
	RandomVariable[] cacheValuesOfUnderlying;

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
			double fixingDate = bermudanOption.getFixingDates()[period];
			double exerciseDate = fixingDate;
			double periodLength = bermudanOption.getPeriodLengths()[period];
			double paymentDate = bermudanOption.getPaymentDates()[period];
			double notional = bermudanOption.getPeriodNotionals()[period];
			double swaprate = bermudanOption.getSwaprates()[period];

			// Get random variables - note that this is the rate at simulation time =
			// exerciseDate
			RandomVariable libor = model.getLIBOR(fixingDate, fixingDate, fixingDate + periodLength);

			// foreach(path) values[path] += notional * (libor.get(path) - swaprate) *
			// periodLength / numeraire.get(path) * monteCarloProbabilities.get(path);
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

			optionValue = triggerValues.choose(continuationValue, exerciseValue);

			exerciseTime = triggerValues.choose(exerciseTime, new Scalar(exerciseDate));
			cacheOptionValues[period] = optionValue;

		}

		// Note that values is a relative price - no numeraire division is required
		RandomVariable numeraireAtZero = model.getNumeraire(evaluationTime);
		RandomVariable monteCarloProbabilitiesAtZero = model.getMonteCarloWeights(evaluationTime);

		optionValue = optionValue.mult(numeraireAtZero).div(monteCarloProbabilitiesAtZero);

		return optionValue;
	}

	protected abstract RandomVariable calculateTriggerValues(int period, double fixingDate,
			LIBORModelMonteCarloSimulationModel model) throws CalculationException;

	public abstract RandomVariable[] getBasisFunctions(double evaluationTime, MonteCarloSimulationModel model)
			throws CalculationException;
}