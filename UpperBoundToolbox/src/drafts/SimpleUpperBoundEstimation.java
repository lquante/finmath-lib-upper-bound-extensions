package drafts;

import java.util.ArrayList;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;

public class SimpleUpperBoundEstimation extends AbstractSimpleBoundEstimation {
	AbstractSimpleBoundEstimation lowerBoundMethod;
	SimplestExerciseStrategy exerciseStrategy;

	/**
	 * @param AbstractSimpleBoundEstimation lowerBoundMethod The lower bound method
	 *                                      to be used as a basis for the upper
	 *                                      bound.
	 */
	public SimpleUpperBoundEstimation(AbstractSimpleBoundEstimation lowerBoundMethod) {
		super();
		this.lowerBoundMethod = lowerBoundMethod;
		exerciseStrategy = new SimplestExerciseStrategy();
	}

	@Override
	protected RandomVariable calculateTriggerValues(int period, double fixingDate,
			LIBORModelMonteCarloSimulationModel model) throws CalculationException {
		// calculate lower bound values for dual method:

		this.bermudanOption.setValuationMethod(lowerBoundMethod);
		this.bermudanOption.getValue(model);
		lowerBoundMethod = (SimpleLowerBoundEstimation) this.bermudanOption.getValuationMethod();
		// Calculate martingale as given in 18.13:
		// M(T_0) = U(T_0)
		ArrayList<RandomVariable> martingaleCache = new ArrayList<RandomVariable>();
		RandomVariable martingale = model.getRandomVariableForConstant(0);

		martingaleCache.add(martingale);
		int numberOfForwardPeriods = this.bermudanOption.getFixingDates().length - 1;
		for (int forwardPeriod = 1; forwardPeriod <= numberOfForwardPeriods; forwardPeriod++) {
			RandomVariable u_i = lowerBoundMethod.getCacheOptionValues()[forwardPeriod];
			RandomVariable foresightBiasCorrection = lowerBoundMethod.getCacheConditionalExpectations()[forwardPeriod];
			martingale = martingaleCache.get(forwardPeriod - 1).add(u_i).sub(foresightBiasCorrection);
			martingaleCache.add(martingale);

		}
		// Calculate upper bound backward
		continuationValue = lowerBoundMethod.getCacheConditionalExpectations()[numberOfForwardPeriods]
				.sub(martingaleCache.get(numberOfForwardPeriods));
		RandomVariable triggerValue = null;

		exerciseValue = lowerBoundMethod.getCacheValuesOfUnderlying()[period].sub(martingaleCache.get(period));
		triggerValue = exerciseStrategy.getExerciseIndicators(continuationValue, exerciseValue);

		return triggerValue;

	}
}
