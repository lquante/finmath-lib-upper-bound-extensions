package upperBoundMethods;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.IntStream;

import lowerBoundMethods.AbstractLowerBoundEstimationInputForUpperBound;
import lowerBoundMethods.SimpleLowerBoundEstimation;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariableFromDoubleArray;
import net.finmath.montecarlo.automaticdifferentiation.RandomVariableDifferentiable;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.ConditionalExpectationEstimator;
import net.finmath.stochastic.RandomVariable;

/**
 * Implements the upper bound estimation using AAD delta estimations as proposed
 * by Joshi and Tang 2014
 * 
 * @author Lennart Quante
 *
 */
public class DeltaHedgingUpperBound extends AbstractUpperBoundEstimation {

	LIBORModelMonteCarloSimulationModel model;

	public DeltaHedgingUpperBound(AbstractLowerBoundEstimationInputForUpperBound lowerBoundMethod,
			double weightOfMartingale) {
		super(lowerBoundMethod, weightOfMartingale);

	}

	@Override
	protected ArrayList<RandomVariable> calculateMartingaleApproximation(int evaluationPeriod,
			LIBORModelMonteCarloSimulationModel model, RandomVariable[] cacheUnderlying,
			RandomVariable[] cacheOptionValues, RandomVariable[] triggerValues) throws CalculationException {
		this.model = model;
		double evaluationTime = model.getLiborPeriod(evaluationPeriod);
		double [] fixingDates = this.bermudanSwaption.getFixingDates();
		int numberOfOptionPeriods = fixingDates.length;

		// initialize martingale as lower bound value for period 0 and 1.

		ArrayList<RandomVariable> martingaleCache = new ArrayList<RandomVariable>();
		RandomVariable martingale = cacheOptionValues[0];
		martingaleCache.add(martingale);
		if (this.bermudanSwaption.getFixingDates().length > 1 && evaluationPeriod + 1 < cacheOptionValues.length) {
			martingale = cacheOptionValues[evaluationPeriod + 1];
			martingaleCache.add(martingale);
		}
		// calculate number of LIBOR periods covered by the fixing dates:
		
		double firstFixingDate= fixingDates[0];
		double lastFixingDate= fixingDates[numberOfOptionPeriods-1];
		int firstLIBORIndex = model.getLiborPeriodIndex(firstFixingDate);
		int lastLIBORIndex = model.getLiborPeriodIndex(lastFixingDate);
		
		// approximate remaining martingale components using delta approximation
		Map<Double, RandomVariable> liborMartingale = (deltaMartingaleApproximation(evaluationTime, cacheOptionValues[0], firstLIBORIndex,lastLIBORIndex));
		// filter values for all fixing dates
		for (int fixingDateIndex=0;fixingDateIndex<numberOfOptionPeriods;fixingDateIndex++)
			martingaleCache.add(liborMartingale.get(fixingDates[fixingDateIndex]));
		// return array of martingale approximations:
		return martingaleCache;
	}

	/**
	 * Method to approximate the martingale process based on delta approximations
	 * from AAD.
	 * 
	 * @param evaluationTime      The time at which the martingale should be
	 *                            evaluated (i.e. discounted to)
	 * @param firstLIBORIndex 		The first LIBOR period for which the martingale shall be approximated.
	 * @param lastLIBORIndex 		The last LIBOR period for which the martingale shall be approximated.
	 * @return An array list with all estimated martingale RandomVariables
	 * @throws CalculationException
	 */
	private Map<Double, RandomVariable> deltaMartingaleApproximation(double evaluationTime,
			RandomVariable cacheOptionValue, int firstLIBORIndex, int lastLIBORIndex) throws CalculationException {

		// Going forward in time we monitor the hedge deltas on each path.
		RandomVariableDifferentiable value = (RandomVariableDifferentiable) cacheOptionValue;
		// Gradient of option value to replicate
		Map<Long, RandomVariable> gradient = value.getGradient();
		
		// initialize martingale cache
		Map<Double, RandomVariable> martingaleCache = new HashMap<>();
		// loop over all discretization dates of the libor discretization to calculate all martingale
		// components
		IntStream.range(firstLIBORIndex, lastLIBORIndex).parallel().forEach(liborTimeIndex ->

		{

			// get current model time index and time
			double liborTime = model.getLiborPeriod(liborTimeIndex);
			// calculate deltas for every fixing date
			RandomVariable[] deltas = null;
			try {
				deltas = getDeltas(gradient, liborTime,firstLIBORIndex,lastLIBORIndex);
			} catch (CalculationException e1) {

				e1.printStackTrace();
			}
			// approximate martingale according to the deltas
			
			RandomVariable martingale = model.getRandomVariableForConstant(0);
			RandomVariable [] martingaleComponents = new RandomVariable[lastLIBORIndex-firstLIBORIndex];
			RandomVariable[] deltaInput = deltas;
			IntStream.range(firstLIBORIndex, lastLIBORIndex).parallel().forEach(liborPeriodIndex->
			//for (int liborPeriodIndex = firstLIBORIndex; liborPeriodIndex < lastLIBORIndex; liborPeriodIndex++) 
					{
				// get times for current rate calculation

				double fixingDate = model.getLiborPeriod(liborPeriodIndex);
				double paymenDate = model.getLiborPeriod(liborPeriodIndex+1);
				double currentPeriodLength = paymenDate-fixingDate;
				// get times for forward calculation
				double forwardFixingDate = paymenDate;
				double forwardPaymenDate = model.getLiborPeriod(liborPeriodIndex+2);
				double forwardPeriodLength = forwardPaymenDate-forwardFixingDate;

				RandomVariable currentRate = null;
				try {
					currentRate = model.getLIBOR(liborTime, fixingDate, paymenDate);
				} catch (CalculationException e) {

					e.printStackTrace();
				}
				RandomVariable bondValue = currentRate.mult(currentPeriodLength);
				RandomVariable forwardAtTimeIndex = null;
				try {
					forwardAtTimeIndex = model.getLIBOR(liborTime, forwardFixingDate, forwardPaymenDate);
				} catch (CalculationException e) {

					e.printStackTrace();
				}
				RandomVariable forwardValue = forwardAtTimeIndex.mult(forwardPeriodLength);

				// calculate martingale according to 5.1 of Joshi / Tang (2014)
				martingaleComponents[liborPeriodIndex-firstLIBORIndex] = (deltaInput[liborPeriodIndex-firstLIBORIndex].mult(forwardValue.sub(bondValue)));
				
			});
			for (int i=0; i<(lastLIBORIndex-firstLIBORIndex);i++)
				martingale.add(martingaleComponents[i]);
			martingaleCache.put(liborTime,martingale);
		});
		return martingaleCache;

	}

	/**
	 * Method to estimate the relevant deltas for the martingale approximation.
	 * 
	 * @param gradient       The already calculated gradient map.
	 * @param modelTimeIndex The model time index for which the deltas should be
	 *                       calculated
	 * @return The delta approximation as a RandomVariable.
	 * @throws CalculationException
	 */
	private RandomVariable[] getDeltas(Map<Long, RandomVariable> gradient, double liborTime, int firstLIBORIndex, int lastLIBORIndex)
			throws CalculationException {

		int numberOfLIBORDates = lastLIBORIndex-firstLIBORIndex;
		RandomVariable[] deltas = new RandomVariable[numberOfLIBORDates];

		// loop parallelized over all fixing dates of the option to calculate each
		// individual delta
		
		IntStream.range(firstLIBORIndex, lastLIBORIndex).parallel().forEach(liborPeriodIndex -> {
			
			SimpleLowerBoundEstimation valuationMethod = (SimpleLowerBoundEstimation) this.bermudanSwaption
					.getValuationMethod();
			RandomVariable exerciseTime = valuationMethod.getExerciseTime();
			Map<Double, Long> liborIDs = valuationMethod.getLiborIDs();
			// get times for current rate calculation
			double fixingDate = model.getLiborPeriod(liborPeriodIndex);
			double paymenDate = model.getLiborPeriod(liborPeriodIndex+1);
			RandomVariable currentRate = null;
			try {
				currentRate = (model.getLIBOR(liborTime, fixingDate, paymenDate));
			} catch (CalculationException e) {

				e.printStackTrace();
			}

			RandomVariable currentNumeraire = null;
			try {
				currentNumeraire = model.getNumeraire(liborTime);
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
				indicator = exerciseTime.sub(liborTime + 0.001).choose(new RandomVariableFromDoubleArray(1.0),
						new RandomVariableFromDoubleArray(0.0));
			}

			// Create a conditional expectation estimator with some basis functions
			// (predictor variables) for conditional expectation estimation.
			ArrayList<RandomVariable> basisFunctions = getRegressionBasisFunctionsBinning(currentRate, indicator);
			ConditionalExpectationEstimator conditionalExpectationOperator = new MonteCarloConditionalExpectationRegression(
					basisFunctions.toArray(new RandomVariable[0]));

			delta = delta.getConditionalExpectation(conditionalExpectationOperator);

			// store delta in delta array
			deltas[liborPeriodIndex-firstLIBORIndex] = delta;
		});
		return deltas;

	}

	/**
	 * Method to provide basis functions for delta approximation using binning.
	 * 
	 * @param underlying        The underlying value of the deltas.
	 * @param exerciseIndicator The exercise indicator of the option.
	 * @return The basis functions to estimate the conditional expectation of the
	 *         deltas.
	 */
	private ArrayList<RandomVariable> getRegressionBasisFunctionsBinning(RandomVariable underlying,
			RandomVariable exerciseIndicator) {
		ArrayList<RandomVariable> basisFunctions = new ArrayList<RandomVariable>();

		if (underlying.isDeterministic()) {
			basisFunctions.add(underlying);
		} else {
			int numberOfBins =  20; //this.bermudanSwaption.getFixingDates().length * 10; optimal number of bins?
			double[] values = underlying.getRealizations();
			Arrays.sort(values);
			IntStream.range(0, numberOfBins).parallel().forEach(i->{	
				double binLeft = values[(int) (((double) i / (double) numberOfBins) * values.length)];
				RandomVariable basisFunction = underlying.sub(binLeft)
						.choose(new RandomVariableFromDoubleArray(1.0), new RandomVariableFromDoubleArray(0.0))
						.mult(exerciseIndicator);
				basisFunctions.add(basisFunction);
			});
		}

		return basisFunctions;
	}
}