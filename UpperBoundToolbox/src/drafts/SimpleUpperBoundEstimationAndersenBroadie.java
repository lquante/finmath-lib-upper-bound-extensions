package drafts;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.function.DoubleBinaryOperator;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;

public class SimpleUpperBoundEstimationAndersenBroadie extends AbstractUpperBoundEstimation {
	AbstractSimpleBoundEstimation lowerBoundMethod;
	SimplestExerciseStrategy exerciseStrategy;
	private int numberOfSubsimulationsStepA;
	private int numberOfSubsimulationsStepB;

	/**
	 * @param AbstractSimpleBoundEstimation lowerBoundMethod The lower bound method
	 *                                      to be used as a basis for the upper
	 *                                      bound.
	 */
	public SimpleUpperBoundEstimationAndersenBroadie(AbstractSimpleBoundEstimation lowerBoundMethod,
			int numberOfSubsimulationsStepA, int numberOfSubsimulationsStepB) {
		super();
		this.lowerBoundMethod = lowerBoundMethod;
		this.numberOfSubsimulationsStepA = numberOfSubsimulationsStepA;
		this.numberOfSubsimulationsStepB = numberOfSubsimulationsStepB;
		exerciseStrategy = new SimplestExerciseStrategy();
	}

	@Override
	protected RandomVariable calculateOptionValue(int period, LIBORModelMonteCarloSimulationModel model)
			throws CalculationException {
		// calculate lower bound values for dual method:
		this.bermudanOption.setValuationMethod(lowerBoundMethod);
		this.bermudanOption.getValue(model);
		lowerBoundMethod = (SimpleLowerBoundEstimation) this.bermudanOption.getValuationMethod();
		RandomVariable[] cacheUnderlying = lowerBoundMethod.getCacheValuesOfUnderlying();
		RandomVariable[] cacheOptionValues = lowerBoundMethod.getCacheOptionValues();

		// Check exercise condition of lower bound method
		RandomVariable triggerValues = lowerBoundMethod.getTriggerValues();
		// create models for subsimulations
		Map<String, Object> modifiedDataA = new HashMap<String, Object>();
		modifiedDataA.put("numberOfPaths", numberOfSubsimulationsStepA);
		Map<String, Object> modifiedDataB = new HashMap<String, Object>();
		modifiedDataB.put("numberOfPaths", numberOfSubsimulationsStepB);
		LIBORModelMonteCarloSimulationModel modelStepA = (LIBORModelMonteCarloSimulationModel) model
				.getCloneWithModifiedData(modifiedDataA);
		LIBORModelMonteCarloSimulationModel modelStepB = (LIBORModelMonteCarloSimulationModel) model
				.getCloneWithModifiedData(modifiedDataB);
		ArrayList<RandomVariable> martingaleCache = new ArrayList<RandomVariable>();
		// initialize martingale as lower bound value for period 0 and 1.
		RandomVariable martingale = cacheOptionValues[0];
		martingaleCache.add(martingale);
		if (this.bermudanOption.getFixingDates().length > 1) {
			martingale = cacheOptionValues[1];
			martingaleCache.add(martingale);
		}

		// determine number of additional martingale components to be estimated
		int numberOfForwardPeriods = this.bermudanOption.getFixingDates().length - 1;
		for (int forwardPeriod = 2; forwardPeriod <= numberOfForwardPeriods; forwardPeriod++) {
			// initialize random variables and models for subsimulation
			RandomVariable subsimulationValue;
			RandomVariable subsimulatedOptionValue;
			BermudanSwaption restartedBermudanA = this.bermudanOption
					.getCloneWithModifiedStartingPeriod(forwardPeriod / 2);
			BermudanSwaption restartedBermudanB = this.bermudanOption
					.getCloneWithModifiedStartingPeriod(forwardPeriod / 2);
			RandomVariable discountFactor = modelStepA.getNumeraire(forwardPeriod / 2)
					.div(modelStepA.getMonteCarloWeights(forwardPeriod));
			// case 2a / b subsimulation conducted via trigger values from lower bound:
			if (forwardPeriod < numberOfForwardPeriods) {
				subsimulationValue = triggerValues.choose(
						restartedBermudanA.getValue(forwardPeriod / 2, modelStepA).div(discountFactor),
						cacheUnderlying[forwardPeriod]);
			} else
				subsimulationValue = triggerValues.choose(modelStepA.getRandomVariableForConstant(0),
						modelStepB.getRandomVariableForConstant(0));
			// check if forwardPeriod<final period
			if (forwardPeriod + 1 < numberOfForwardPeriods) {
				RandomVariable valueCorrectionTerm = restartedBermudanB.getValue(forwardPeriod / 2, modelStepB)
						.div(discountFactor).sub(cacheUnderlying[forwardPeriod]);
				subsimulatedOptionValue = triggerValues.choose(modelStepA.getRandomVariableForConstant(0),
						valueCorrectionTerm);

			} else
				subsimulatedOptionValue = triggerValues.choose(modelStepA.getRandomVariableForConstant(0),
						modelStepB.getRandomVariableForConstant(0));

			// adjust martingale and save in martingale cache
			martingale = martingaleCache.get(forwardPeriod - 1).add(subsimulationValue)
					.sub(cacheOptionValues[forwardPeriod - 1]).sub(subsimulatedOptionValue);
			martingaleCache.add(martingale);

		}
		// final upper bound value calculated forward via martingale (not like in simple
		// upper bound, no trigger values)
		RandomVariable upperBoundValue = cacheUnderlying[0].sub(martingaleCache.get(0));
		// calculate the maximum from the simulated martingale for each remaining period
		for (int forwardPeriod = 1; forwardPeriod <= numberOfForwardPeriods; forwardPeriod++) {
			DoubleBinaryOperator max = (x, y) -> Math.max(x, y);
			upperBoundValue = upperBoundValue.apply(max,
					cacheUnderlying[forwardPeriod].sub(martingaleCache.get(forwardPeriod)));
		}
		// discount and adjust for monte carlo weights, adjust by martingale estimation
		// (Delta in Andersen Broadie paper)
		upperBoundValue = cacheOptionValues[0].mult(model.getNumeraire(0)).div(model.getMonteCarloWeights(0))
				.add(upperBoundValue.getAverage());
		return upperBoundValue;
	}
}
