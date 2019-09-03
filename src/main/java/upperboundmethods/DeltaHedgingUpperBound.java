package upperboundmethods;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.IntStream;

import lowerboundmethods.AbstractLowerBoundEstimationInputForUpperBound;
import lowerboundmethods.SimpleLowerBoundEstimation;
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

	/**
	 * @param lowerBoundMethod         The lower bound method to be used as a basis
	 *                                 for the upper bound.
	* @param weightOfMartingale        The weighting scheme for point value approximation - 0 = lower bound, 1= upper bound, 0.5 = A-B point wise estimate.	    
	 */
	public DeltaHedgingUpperBound(AbstractLowerBoundEstimationInputForUpperBound lowerBoundMethod,
			double weightOfMartingale) {
		super(lowerBoundMethod, weightOfMartingale);

	}
	/**Constructs DeltaHedging upper bound estimator with fixed simple lower bound and martingale weight of 1.
	 * 
	 */
	public DeltaHedgingUpperBound() {
		super(new SimpleLowerBoundEstimation(),1);
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

		double firstFixingDate= fixingDates[2];
		double lastFixingDate= fixingDates[numberOfOptionPeriods-1];
		int firstLIBORIndex = model.getLiborPeriodIndex(firstFixingDate);
		int lastLIBORIndex = model.getLiborPeriodIndex(lastFixingDate);

		// approximate remaining martingale components using delta approximation
		Map<Double, RandomVariable> liborMartingale = (deltaMartingaleApproximation(evaluationTime, cacheOptionValues[0], firstLIBORIndex,lastLIBORIndex));
		// filter values for all fixing dates
		for (int fixingDateIndex=2;fixingDateIndex<numberOfOptionPeriods;fixingDateIndex++)
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
		IntStream.range(firstLIBORIndex, lastLIBORIndex+1).parallel().forEach(liborTimeIndex ->
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
			{
				RandomVariable bondValue = calculateOnePeriodBondValue(liborTime,liborPeriodIndex); 
				RandomVariable forwardValue = calculateOnePeriodBondValue(liborTime,liborPeriodIndex+1); 
				// calculate martingale according to 5.1 of Joshi / Tang (2014)
				martingaleComponents[liborPeriodIndex-firstLIBORIndex]= (deltaInput[liborPeriodIndex-firstLIBORIndex].mult(forwardValue.sub(bondValue)));
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
		RandomVariable[] deltas = new RandomVariable[lastLIBORIndex-firstLIBORIndex];
		// loop parallelized over all fixing dates of the option to calculate each individual delta
		IntStream.range(firstLIBORIndex, lastLIBORIndex).parallel().forEach(liborPeriodIndex -> {
			SimpleLowerBoundEstimation valuationMethod = (SimpleLowerBoundEstimation) this.bermudanSwaption
					.getValuationMethod();
			// get fixing date from LIBOR index 
			double fixingDate = model.getLiborPeriod(liborPeriodIndex);
			// get current LIBOR and current numeraire
			RandomVariable currentRate = null;
			RandomVariable currentNumeraire = null;
			try {
				currentNumeraire = model.getNumeraire(fixingDate);
				currentRate = model.getLIBOR(model.getTimeIndex(liborTime), liborPeriodIndex);
			} catch (CalculationException e1) {
				e1.printStackTrace();
			}
			// get gradient with respect to the forward rate and adjust by numeraire
			RandomVariable delta = gradient.get(valuationMethod.getLiborIDs().get(fixingDate));
			
			if (delta == null) {
				delta = currentRate.mult(0.0);
			}
			// adjust by numeraire
			delta= delta.mult(currentNumeraire);
			// get exerciseIndicator
			RandomVariable indicator = calculateExerciseIndicator  (liborTime,valuationMethod.getExerciseTime());
			// Create a conditional expectation estimator with some basis functions
			// (predictor variables) for conditional expectation estimation.
			ArrayList<RandomVariable> basisFunctions = getRegressionBasisFunctionsBinning(currentRate, indicator);
			ConditionalExpectationEstimator conditionalExpectationOperator = new MonteCarloConditionalExpectationRegression(
					basisFunctions.toArray(new RandomVariable[0]));
			// calculate conditional expectation and store delta in delta array
			deltas[liborPeriodIndex-firstLIBORIndex] = delta.getConditionalExpectation(conditionalExpectationOperator);
		});
		return deltas;
	}

	/**
	 * Method to calculate exercise indicator for a given time and the exerciseTime from a valuation
	 * @param time the time at which the exercise indicator should be evaluated
	 * @param exerciseTime the exercise time from the valuationMethod of the Bermudan
	 * @return the exercise indicator as random variable
	 */
	private RandomVariable calculateExerciseIndicator(double time, RandomVariable exerciseTime) {
		RandomVariable indicator = new RandomVariableFromDoubleArray(1.0);
		if (exerciseTime != null) {
			indicator = exerciseTime.sub(time + 0.001).choose(new RandomVariableFromDoubleArray(1.0),
					new RandomVariableFromDoubleArray(0.0));
		}
		return indicator;
	}
	
	/** 
	 * Small calculation method for the value of a simple bond
	 * @param evaluationTime time to evaluate the primitive bond
	 * @param liborPeriodIndex The LIBOR index of the rate to be used
	 * @return the value of a bond with notional 1 over 1 period of the LIBOR model
	 */
	private RandomVariable calculateOnePeriodBondValue(double evaluationTime, int liborPeriodIndex) {
		// get times for rate calculation
		double fixingDate = model.getLiborPeriod(liborPeriodIndex);
		double paymenDate = model.getLiborPeriod(liborPeriodIndex+1);
		double periodLength = paymenDate-fixingDate;
		// get LIBOR		
		RandomVariable rate = null;
		try {
			rate = model.getLIBOR(evaluationTime, fixingDate, paymenDate);
		} catch (CalculationException e) {

			e.printStackTrace();
		}
		// return one period bond value		
		return rate.mult(periodLength);
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
			int numberOfBins =  20; //optimal number of bins?
			double[] values = underlying.getRealizations();
			Arrays.sort(values);
			for (int i =0; i<numberOfBins;i++){	
				double binLeft = values[(int) (((double) i / (double) numberOfBins) * values.length)];
				RandomVariable basisFunction = underlying.sub(binLeft)
						.choose(new RandomVariableFromDoubleArray(1.0), new RandomVariableFromDoubleArray(0.0))
						.mult(exerciseIndicator);
				basisFunctions.add(basisFunction);
			}
		}
		return basisFunctions;
	}
}