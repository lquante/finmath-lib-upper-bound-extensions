package upperBoundMethods;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.stream.IntStream;

import lowerBoundMethods.AbstractLowerBoundEstimationInputForUpperBound;
import lowerBoundMethods.AbstractLowerBoundEstimationWithoutCaching;
import lowerBoundMethods.SimpleLowerBoundEstimation;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariableFromDoubleArray;
import net.finmath.montecarlo.automaticdifferentiation.RandomVariableDifferentiable;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.ConditionalExpectationEstimator;
import net.finmath.stochastic.RandomVariable;

/**
 * Implements the upper bound estimation using AAD delta estimations as proposed by Joshi and Tang 2014
 * @author Lennart Quante
 *
 */
public class DeltaHedgingUpperBound extends AbstractUpperBoundEstimation {

	LIBORModelMonteCarloSimulationModel model;

	public DeltaHedgingUpperBound(AbstractLowerBoundEstimationInputForUpperBound lowerBoundMethod) {
		super(lowerBoundMethod);

	}

	@Override
	protected double calculateMartingaleApproximation(int evaluationPeriod, LIBORModelMonteCarloSimulationModel model,
			RandomVariable[] cacheUnderlying, RandomVariable[] cacheOptionValues, RandomVariable[] triggerValues)
					throws CalculationException {
		this.model = model;
		double evaluationTime =model.getTime(evaluationPeriod);
		int numberOfOptionPeriods = this.bermudanOption.getFixingDates().length;

		// initialize martingale as lower bound value for period 0 and 1.

		ArrayList<RandomVariable> martingaleCache = new ArrayList<RandomVariable>();
		RandomVariable martingale = cacheOptionValues[0];
		martingaleCache.add(martingale);
		if (this.bermudanOption.getFixingDates().length > 1 && evaluationPeriod + 1 < cacheOptionValues.length) {
			martingale = cacheOptionValues[evaluationPeriod + 1];
			martingaleCache.add(martingale);
		}
		// approximate remaining martingale components using delta approximation

		martingaleCache.addAll(deltaMartingaleApproximation(evaluationTime, numberOfOptionPeriods));

		// calculate the maximum from the simulated martingale process for each remaining period

		RandomVariable martingaleApproximation = model.getRandomVariableForConstant(0);
		for (int forwardPeriod = 0; forwardPeriod < martingaleCache.size()
				&& forwardPeriod < numberOfOptionPeriods; forwardPeriod++)
			martingaleApproximation = martingaleApproximation.floor(cacheUnderlying[forwardPeriod].sub(martingaleCache.get(forwardPeriod)));

		// return average of martingale process
		return  martingaleApproximation.getAverage();

	}


	/**
	 * Method to approximate the martingale process based on delta approximations from AAD.
	 * @param evaluationTime The time at which the martingale should be evaluated (i.e. discounted to)
	 * @param numberOfFixingDates The number of periods for which the martingale process needs to be approximated.
	 * @return An array list with all estimated martingale RandomVariables
	 * @throws CalculationException
	 */
	private ArrayList<RandomVariable> deltaMartingaleApproximation(double evaluationTime, int numberOfFixingDates)
			throws CalculationException {

		/*
		 * Going forward in time we monitor the hedge deltas on each path.
		 */

		long timingValuationStart = System.currentTimeMillis();

		RandomVariableDifferentiable value = (RandomVariableDifferentiable) this.bermudanOption.getValue(evaluationTime, model);

		long timingValuationEnd = System.currentTimeMillis();

		// Gradient of option value to replicate
		long timingDerivativeStart = System.currentTimeMillis();
		Map<Long, RandomVariable> gradient = value.getGradient();
		long timingDerivativeEnd = System.currentTimeMillis();


		// loop over all discretization dates of the option to calculate all martingale
		// components
		// initialize martingale cache
		ArrayList<RandomVariable> martingaleCache = new ArrayList<RandomVariable>();

		IntStream.range(2,numberOfFixingDates).parallel().forEach(optionTimeIndex->


		{

			// get current model time index and time
			int modelTimeIndex = model.getTimeIndex(this.bermudanOption .getFixingDates()[optionTimeIndex]);
			double modelTime= model.getTime(modelTimeIndex);
			// calculate deltas for every fixing date
			RandomVariable[] deltas = null;
			try {
				deltas = getDeltas(gradient, modelTime);
			} catch (CalculationException e1) {

				e1.printStackTrace();
			}
			/*
			 * calculate the martingale according to the deltas
			 */
			RandomVariable martingale = model.getRandomVariableForConstant(0);
			
			
			for (int fixingDateIndex = 0; fixingDateIndex < numberOfFixingDates-1; fixingDateIndex++) {
				// get times for current rate calculation
				double fixingDate = this.bermudanOption.getFixingDates()[fixingDateIndex];
				double paymenDate = this.bermudanOption.getPaymentDates()[fixingDateIndex];
				double currentPeriodLength = this.bermudanOption.getPeriodLengths()[fixingDateIndex];
				// get times for forward calculation
				double forwardFixingDate = this.bermudanOption.getFixingDates()[fixingDateIndex+1];
				double forwardPaymenDate = this.bermudanOption.getPaymentDates()[fixingDateIndex+1];
				double forwardPeriodLength = this.bermudanOption.getPeriodLengths()[fixingDateIndex+1]; 

				RandomVariable currentRate = null;
				try {
					currentRate = model.getLIBOR(modelTime, fixingDate,
							paymenDate);
				} catch (CalculationException e) {

					e.printStackTrace();
				}
				RandomVariable bondValue = currentRate.mult(currentPeriodLength);
				RandomVariable forwardAtTimeIndex = null;
				try {
					forwardAtTimeIndex = model.getLIBOR(modelTime,forwardFixingDate ,
							forwardPaymenDate);
				} catch (CalculationException e) {

					e.printStackTrace();
				}
				RandomVariable forwardValue = forwardAtTimeIndex.mult(forwardPeriodLength);

				// calculate martingale according to 5.1 of Joshi / Tang (2014)
				martingale = martingale.add(deltas[fixingDateIndex].mult(forwardValue.sub(bondValue)));
			}
			martingaleCache.add(martingale);
		}
				);
		return martingaleCache;

	}

	/**
	 * Method to estimate the relevant deltas for the martingale approximation.
	 * @param gradient The already calculated gradient map. 
	 * @param modelTimeIndex The model time index for which the deltas should be calculated
	 * @return The delta approximation as a RandomVariable.
	 * @throws CalculationException
	 */
	private RandomVariable[] getDeltas(Map<Long, RandomVariable> gradient, double modelTime) throws CalculationException {

		int numberOfFixingDates = this.bermudanOption.getFixingDates().length;
		RandomVariable[] deltas = new RandomVariable[numberOfFixingDates];

		// loop parallelized over all fixing dates of the option to calculate each individual delta

		IntStream.range(0,numberOfFixingDates).parallel().forEach(i-> {

			SimpleLowerBoundEstimation valuationMethod = (SimpleLowerBoundEstimation) this.bermudanOption
					.getValuationMethod();
			RandomVariable exerciseTime = valuationMethod.getExerciseTime();
			Map<Double, Long> liborIDs = valuationMethod.getLiborIDs();
			// get times for current rate calculation
			double fixingDate = this.bermudanOption.getFixingDates()[i];
			double paymenDate = this.bermudanOption.getPaymentDates()[i];

			RandomVariable currentRate = null;
			try {
				currentRate = ((RandomVariableDifferentiable) model.getLIBOR(modelTime, fixingDate,
						paymenDate));
			} catch (CalculationException e) {
				
				e.printStackTrace();
			}

			RandomVariable currentNumeraire = null;
			try {
				currentNumeraire = model.getNumeraire(modelTime);
			} catch (CalculationException e) {
				
				e.printStackTrace();
			}
			// get gradient with respect to the forward rate (liborID's are calculated in
			// backward algorithm thus inversion of index)
			RandomVariable delta = gradient.get(liborIDs.get(fixingDate));

			if (delta == null) {
				delta = currentRate.mult(0.0);
			}

			// adjust by numeraire
			delta = delta.mult(currentNumeraire);
			// get exerciseIndicator
			RandomVariable indicator = new RandomVariableFromDoubleArray(1.0);
			if (exerciseTime != null) {
				indicator = exerciseTime.sub(modelTime + 0.001)
						.choose(new RandomVariableFromDoubleArray(1.0), new RandomVariableFromDoubleArray(0.0));
			}

			// Create a conditional expectation estimator with some basis functions
			// (predictor variables) for conditional expectation estimation.
			ArrayList<RandomVariable> basisFunctions = getRegressionBasisFunctionsBinning(currentRate,
					indicator);
			ConditionalExpectationEstimator conditionalExpectationOperator = new MonteCarloConditionalExpectationRegression(
					basisFunctions.toArray(new RandomVariable[0]));

			delta = delta.getConditionalExpectation(conditionalExpectationOperator);

			// store delta in delta array
			deltas[i] = delta;
		}
				);
		return deltas;

	}



	/**
	 * Method to provide basis functions for delta approximation using binning.
	 * @param underlying The underlying value of the deltas.
	 * @param exerciseIndicator The exercise indicator of the option.
	 * @return The basis functions to estimate the conditional expectation of the deltas.
	 */
	private ArrayList<RandomVariable> getRegressionBasisFunctionsBinning(RandomVariable underlying,
			RandomVariable exerciseIndicator) {
		ArrayList<RandomVariable> basisFunctions = new ArrayList<RandomVariable>();

		if (underlying.isDeterministic()) {
			basisFunctions.add(underlying);
		} else {
			int numberOfBins = 20;
			double[] values = underlying.getRealizations();
			Arrays.sort(values);
			IntStream.range( 0,  numberOfBins).parallel().forEach(i-> {
				double binLeft = values[(int) ((double) i / (double) numberOfBins * values.length)];
				RandomVariable basisFunction = underlying.sub(binLeft)
						.choose(new RandomVariableFromDoubleArray(1.0), new RandomVariableFromDoubleArray(0.0))
						.mult(exerciseIndicator);
				basisFunctions.add(basisFunction);
			}
					);
		}

		return basisFunctions;
	}
}