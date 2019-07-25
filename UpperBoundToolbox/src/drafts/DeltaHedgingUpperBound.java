package drafts;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Set;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.products.Forward;
import net.finmath.montecarlo.RandomVariableFromDoubleArray;
import net.finmath.montecarlo.assetderivativevaluation.products.BermudanOption;
import net.finmath.montecarlo.automaticdifferentiation.RandomVariableDifferentiable;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.ConditionalExpectationEstimator;
import net.finmath.stochastic.RandomVariable;

public class DeltaHedgingUpperBound extends AbstractUpperBoundEstimation {

	LIBORModelMonteCarloSimulationModel model;

	public DeltaHedgingUpperBound(AbstractSimpleBoundEstimation lowerBoundMethod) {
		super(lowerBoundMethod);
		
	}

	@Override
	protected RandomVariable calculateOptionValue(int period, LIBORModelMonteCarloSimulationModel model,
			RandomVariable[] cacheUnderlying, RandomVariable[] cacheOptionValues, RandomVariable[] triggerValues)
			throws CalculationException {
		this.model = model;
		double evaluationTime = 0;
		return deltaApproximation(evaluationTime, cacheUnderlying, cacheOptionValues);
	}

	// 1st step: delta hedge approximation method

	private RandomVariable deltaApproximation(double evaluationTime, RandomVariable[] cacheUnderlying,
			RandomVariable[] cacheOptionValues) throws CalculationException {
		// Ask the model for its discretization
		int timeIndexEvaluationTime = model.getTimeIndex(evaluationTime);

		/*
		 * Going forward in time we monitor the hedge portfolio on each path.
		 */

		long timingValuationStart = System.currentTimeMillis();

		RandomVariableDifferentiable value = (RandomVariableDifferentiable) this.bermudanOption.getValue(0, model);

		long timingValuationEnd = System.currentTimeMillis();

		RandomVariable valueOfOption = model.getRandomVariableForConstant(value.getAverage());

		// Initialize the portfolio to zero bonds and as much cash as the lower bound
		// valuation predicts.

		RandomVariable numeraireToday = model.getNumeraire(0.0);

		// We store the composition of the hedge portfolio (depending on the path)
		RandomVariable amountOfNumeraireAsset = valueOfOption.div(numeraireToday);
		int numberOfPeriods = this.bermudanOption.getFixingDates().length;
		RandomVariable[] amountOfUnderlyingForwards = new RandomVariable[numberOfPeriods];
		for (int i = 0; i < numberOfPeriods; i++)
			amountOfUnderlyingForwards[i] = model.getRandomVariableForConstant(0.0);

		// Gradient of option value to replicate
		long timingDerivativeStart = System.currentTimeMillis();
		Map<Long, RandomVariable> gradient = value.getGradient();
		long timingDerivativeEnd = System.currentTimeMillis();

		double lastOperationTimingValuation = (timingValuationEnd - timingValuationStart) / 1000.0;
		double lastOperationTimingDerivative = (timingDerivativeEnd - timingDerivativeStart) / 1000.0;
		
		// calculate deltas for every time point between current time (0) and evaluationTime
		for (int timeIndex = 0; timeIndex < timeIndexEvaluationTime; timeIndex++) {

			RandomVariable[] deltas = getDeltas(gradient, timeIndex);

			/*
			 * Change the portfolio according to the trading strategy
			 */
			// loop over all fixing dates of the option to calculate each individual hedge

			for (int i = 0; i < numberOfPeriods; i++) {
				/*
				 * Change the portfolio according to the trading strategy
				 */
				// get times for forward calculation
				double fixingDate = this.bermudanOption.getFixingDates()[i];
				double periodLength = this.bermudanOption.getPeriodLengths()[i];

				RandomVariable forwardAtTimeIndex = model.getLIBOR(model.getTime(timeIndex), fixingDate,
						fixingDate + periodLength);
				
				RandomVariable forwardValue = forwardAtTimeIndex.mult(periodLength);
				RandomVariable numeraireAtTimeIndex = model.getNumeraire(timeIndex);

				// Determine the delta hedge
				RandomVariable newNumberOfForwards = deltas[i];
				RandomVariable bondsToBuy = newNumberOfForwards.sub(amountOfUnderlyingForwards[i]);

				// Ensure self financing
				RandomVariable numeraireAssetsToSell = bondsToBuy.mult(forwardValue).div(numeraireAtTimeIndex);
				RandomVariable newNumberOfNumeraireAsset = amountOfNumeraireAsset.sub(numeraireAssetsToSell);

				// Update portfolio
				amountOfNumeraireAsset = newNumberOfNumeraireAsset;

				amountOfUnderlyingForwards[i] = newNumberOfForwards;
			}

		}

		/*
		 * At maturity, calculate the value of the replication portfolio
		 */

		// Get value of underlying and numeraire assets
		RandomVariable numeraireAtEvaluationTime = model.getNumeraire(evaluationTime);
		RandomVariable portfolioValue = amountOfNumeraireAsset.mult(numeraireAtEvaluationTime);

		for (int i = 0; i < numberOfPeriods; i++) {
			double fixingDate = this.bermudanOption.getFixingDates()[i];
			double periodLength = this.bermudanOption.getPeriodLengths()[i];

			RandomVariable forwardValueAtEvaluationTime = ((RandomVariableDifferentiable) model.getLIBOR(evaluationTime,
					fixingDate, fixingDate + periodLength)).mult(periodLength);
			portfolioValue = portfolioValue.add(amountOfUnderlyingForwards[i].mult(forwardValueAtEvaluationTime));
		}

		return portfolioValue;

	}

	private RandomVariable[] getDeltas(Map<Long, RandomVariable> gradient, int timeIndex) throws CalculationException {

		int numberOfPeriods = this.bermudanOption.getFixingDates().length;
		RandomVariable[] deltas = new RandomVariable[numberOfPeriods];

		// loop over all fixing dates of the option to calculate each individual delta

		for (int i = 0; i < numberOfPeriods; i++) {

			SimpleLowerBoundEstimation valuationMethod = (SimpleLowerBoundEstimation) this.bermudanOption
					.getValuationMethod();
			RandomVariable exerciseTime = valuationMethod.getExerciseTime();
			Set<Long> liborIDs = valuationMethod.getLiborIDs();
			// get times for forward calculation
			double fixingDate = this.bermudanOption.getFixingDates()[i];
			double periodLength = this.bermudanOption.getPeriodLengths()[i];

			RandomVariable forwardAtTimeIndex = ((RandomVariableDifferentiable) model.getLIBOR(model.getTime(timeIndex),
					fixingDate, fixingDate + periodLength));
			RandomVariable numeraireAtTimeIndex = model.getNumeraire(timeIndex);
			Set<Long> test = gradient.keySet();
			// get gradient with respect to the forward rate (liborID's are calculated in backward algorithm thus inversion of index)
			Long forwardID = ((RandomVariableDifferentiable) forwardAtTimeIndex).getID();
			RandomVariable delta = gradient.get(liborIDs.toArray()[numberOfPeriods-1-i]);

			if (delta == null) {
				delta = forwardAtTimeIndex.mult(0.0);
			}

			// adjust by numeraire
			delta = delta.mult(numeraireAtTimeIndex);
			// get exerciseIndicator
			RandomVariable indicator = new RandomVariableFromDoubleArray(1.0);
			if (exerciseTime != null) {
				indicator = exerciseTime.sub(model.getTime(timeIndex) + 0.001)
						.choose(new RandomVariableFromDoubleArray(1.0), new RandomVariableFromDoubleArray(0.0));
			}

			// Create a conditional expectation estimator with some basis functions
			// (predictor variables) for conditional expectation estimation.
			ArrayList<RandomVariable> basisFunctions = getRegressionBasisFunctionsBinning(forwardAtTimeIndex,
					indicator);
			ConditionalExpectationEstimator conditionalExpectationOperator = new MonteCarloConditionalExpectationRegression(
					basisFunctions.toArray(new RandomVariable[0]));

			delta = delta.getConditionalExpectation(conditionalExpectationOperator);

			// store delta in delta array
			deltas[i] = delta;
		}
		return deltas;

	}

	private RandomVariable[] getForwardWeights(Map<Long, RandomVariable> gradient, int timeIndex,
			RandomVariable[] amountOfUnderlyingForwards) throws CalculationException {

		return amountOfUnderlyingForwards;

	}

	private ArrayList<RandomVariable> getRegressionBasisFunctionsBinning(RandomVariable underlying,
			RandomVariable indicator) {
		ArrayList<RandomVariable> basisFunctions = new ArrayList<RandomVariable>();

		if (underlying.isDeterministic()) {
			basisFunctions.add(underlying);
		} else {
			int numberOfBins = 20;
			double[] values = underlying.getRealizations();
			Arrays.sort(values);
			for (int i = 0; i < numberOfBins; i++) {
				double binLeft = values[(int) ((double) i / (double) numberOfBins * values.length)];
				RandomVariable basisFunction = underlying.sub(binLeft)
						.choose(new RandomVariableFromDoubleArray(1.0), new RandomVariableFromDoubleArray(0.0))
						.mult(indicator);
				basisFunctions.add(basisFunction);
			}
		}

		return basisFunctions;
	}

}
