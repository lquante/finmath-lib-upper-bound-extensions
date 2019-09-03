package upperboundmethods;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;

import bermudanswaptionframework.BermudanSwaption;
import bermudanswaptionframework.BermudanSwaptionValueEstimatorInterface;
import bermudanswaptionframework.SimplestExerciseStrategy;
import lowerboundmethods.AbstractLowerBoundEstimationInputForUpperBound;
import lowerboundmethods.AbstractLowerBoundEstimationWithoutCaching;
import lowerboundmethods.SimpleLowerBoundEstimation;
import lowerboundmethods.SimpleLowerBoundEstimationWithoutCaching;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.MonteCarloSimulationModel;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromDoubleArray;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAADFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.LIBORMonteCarloSimulationFromLIBORModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelFromGivenMatrix;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

/**
 * Implements the Andersen-Broadie (2004) method of martingale approximation using
 * sub-simulations.
 * 
 * @author Lennart Quante
 * @version 1.0
 * 
 */
public class AndersenBroadieUpperBoundEstimation extends AbstractUpperBoundEstimation {
	AbstractLowerBoundEstimationWithoutCaching lowerBoundMethod;
	SimplestExerciseStrategy exerciseStrategy;
	private int pathsSubsimulationsStepA;
	private int pathsSubsimulationsStepB;

	/**
	 * @param lowerBoundMethod         The lower bound method to be used as a basis
	 *                                 for the upper bound.
	* @param weightOfMartingale        The weighting scheme for point value approximation - 0 = lower bound, 1= upper bound, 0.5 = A-B point wise estimate.	                        
	 * @param pathsSubsimulationsStepA number of subsimulation Paths in case 2a of
	 *                                 A-B algorithm, i.e. if exercise at current
	 *                                 simulation time.
	 * @param pathsSubsimulationsStepB number of subsimulation Paths in case 2b of
	 *                                 A-B algorithm, i.e. if no exercise at current
	 *                                 simulation time.
	 */

	/**
	 * @param lowerBoundMethod
	 
	 * @param pathsSubsimulationsStepA
	 * @param pathsSubsimulationsStepB
	 */
	public AndersenBroadieUpperBoundEstimation(AbstractLowerBoundEstimationInputForUpperBound lowerBoundMethod, double weightOfMartingale,
			int pathsSubsimulationsStepA, int pathsSubsimulationsStepB) {
		super(lowerBoundMethod,weightOfMartingale);
		this.pathsSubsimulationsStepA = pathsSubsimulationsStepA;
		this.pathsSubsimulationsStepB = pathsSubsimulationsStepB;
		// fixed exercise strategy, maybe add other implementations
		exerciseStrategy = new SimplestExerciseStrategy();
	}
	/**
	 * Creates an AndersenBroadieUpperBound estimator with fixed SimpleLowerBoundEstimation as input, weight of martingale =1 and simplest exercise strategy.
	 * @param subsimulationPaths The number of paths to be used in both subsimulation possibilites
	 */
	public AndersenBroadieUpperBoundEstimation(int subsimulationPaths) {
		super(new SimpleLowerBoundEstimation(),1);
		this.pathsSubsimulationsStepA = subsimulationPaths;
		this.pathsSubsimulationsStepB = subsimulationPaths;
		// fixed exercise strategy, maybe add other implementations
		exerciseStrategy = new SimplestExerciseStrategy();
	}
	
	@Override
	protected ArrayList<RandomVariable> calculateMartingaleApproximation(int evaluationPeriod, LIBORModelMonteCarloSimulationModel model,
			RandomVariable[] cacheUnderlying, RandomVariable[] cacheOptionValues, RandomVariable[] triggerValues)
					throws CalculationException {
		int numberOfPaths = model.getNumberOfPaths();
		// determine number of martingale components to be estimated
		int numberOfOptionPeriods = this.bermudanSwaption.getFixingDates().length;
		// initialize martingale as lower bound value for period 0 and 1.
		ArrayList<RandomVariable> martingaleCache = new ArrayList<RandomVariable>();
		martingaleCache.add(cacheOptionValues[evaluationPeriod]);
		if (evaluationPeriod + 1 < cacheOptionValues.length)
			martingaleCache.add(cacheOptionValues[evaluationPeriod + 1]);
		// estimate martingale for all remaining fixing dates
		for (int martingaleIndex = 2; martingaleIndex < numberOfOptionPeriods; martingaleIndex++) {
			// initialize values for subsimulation
			double currentFixingDate = this.bermudanSwaption.getFixingDates()[martingaleIndex];
			int liborPeriod = model.getLiborPeriodIndex(currentFixingDate);
			// Arrays to store pathwise results to be transformed to RandomVariable
			double[] discountedExerciseValueArray  = new double[numberOfPaths];
			double[] discountedFutureExerciseValueArray = new double[numberOfPaths];
			double[] previousExerciseIndicatorArray = new double[numberOfPaths];
			// Parallelized implementation of pathwise part of A-B algorithm:
			int optionPeriodForInput = martingaleIndex; // to avoid non-final problems in nested environment
			IntStream.range(0, numberOfPaths).parallel().forEach(path -> {
				double[] subsimulationResults = null;
				try {
					subsimulationResults = executeSubsimulation(model, path, optionPeriodForInput, liborPeriod,
							triggerValues, numberOfOptionPeriods, cacheUnderlying);
				} catch (CalculationException e)
				{e.printStackTrace();}
				discountedExerciseValueArray[path] = subsimulationResults[0];
				discountedFutureExerciseValueArray[path] = subsimulationResults[1];
				previousExerciseIndicatorArray[path] = subsimulationResults[2];
			});
			RandomVariable discountedExerciseValue = new RandomVariableFromDoubleArray(currentFixingDate,discountedExerciseValueArray);
			RandomVariable discountedFutureExerciseValue = new RandomVariableFromDoubleArray(currentFixingDate,discountedFutureExerciseValueArray);
			RandomVariable previousExerciseIndicator = new RandomVariableFromDoubleArray(currentFixingDate,previousExerciseIndicatorArray);
			// calculate value of martingale estimation
			martingaleCache.add(martingaleCache.get(martingaleIndex - 1).sub(cacheOptionValues[martingaleIndex - 1])
					.add(discountedExerciseValue)
					.sub((discountedFutureExerciseValue.add(cacheUnderlying[martingaleIndex - 1])).mult(previousExerciseIndicator)));
		}
		return martingaleCache;
	}

	/**
	 * Method to be called by parallelized execution interface for subsimulations
	 * 
	 * @param model
	 * @param path
	 * @param optionPeriod
	 * @param liborPeriod
	 * @param currentFixingDate
	 * @param triggerValues
	 * @param numberOfOptionPeriods
	 * @param cacheUnderlying
	 * @param martingaleIndex
	 * @return The triple
	 *         {discountedExerciseValue,discountedFutureExerciseValue,previousExerciseIndicator}
	 * @throws CalculationException
	 */
	double[] executeSubsimulation(LIBORModelMonteCarloSimulationModel model, int path, int optionPeriod,
			int liborPeriod , RandomVariable[] triggerValues, int numberOfOptionPeriods,
			RandomVariable[] cacheUnderlying) throws CalculationException {
		// initialize variables
		double discountedExerciseValue = 0; double discountedFutureExerciseValue = 0;
		// check trigger values for each path of the main simulation:
		if (triggerValues[optionPeriod].get(path) >= 0) // case 2a
		{
			double currentFixingDate = model.getLiborPeriod(liborPeriod);
			discountedExerciseValue = calculateSubsimulationValue(model,optionPeriod,1000,liborPeriod,1000,path,currentFixingDate, pathsSubsimulationsStepA);
		}
		else // case 2b
		{
			int terminationPeriod = optionPeriod;
			// determine termination period: as soon as triggerValues[path+i] >=0 for i>=optionPeriod+1
			while (terminationPeriod < numberOfOptionPeriods - 1) {
				if (triggerValues[terminationPeriod].get(path) < 0)
					terminationPeriod += 1;
				else
					break;
			}
			// transform to termination time and get LIBORindex
			double terminationTime = this.bermudanSwaption.getFixingDates()[terminationPeriod];
			int terminationTimeIndex = model.getLiborPeriodIndex(terminationTime);
			// retrieve discounted option value from cache:
			discountedExerciseValue = cacheUnderlying[optionPeriod].get(path);
			if (liborPeriod+2 < terminationTimeIndex) {
				// create model for subsimulations 
				double futureExerciseTime = model.getLiborPeriod(liborPeriod+1);
				discountedFutureExerciseValue = calculateSubsimulationValue(model,optionPeriod,terminationPeriod-1,liborPeriod,terminationTimeIndex,path,futureExerciseTime, pathsSubsimulationsStepB);
			}
		}
		// determine previous exercise indicator for martingale calculation
		double previousExerciseIndicator = 0;
		if (triggerValues[optionPeriod - 1].get(path) >= 0)
			previousExerciseIndicator = 1;
		// return results as a triple
		return new double[] { discountedExerciseValue, discountedFutureExerciseValue, previousExerciseIndicator };
	}

	private double calculateSubsimulationValue(LIBORModelMonteCarloSimulationModel model, int startingOptionPeriod,int terminationOptionPeriod, int startingTimeIndex, int terminationTimeIndex,
			int path, double evaluationTime, int subsimulationPaths) throws CalculationException {
		// create model for subsimulation
		LIBORModelMonteCarloSimulationModel modelSubsimulation = createSubsimulationModelTerminating(model, startingTimeIndex,terminationTimeIndex, path,
				subsimulationPaths);
		// create lower bound valuation method without caching, shift option starting date and set this method for the option
		BermudanSwaptionValueEstimatorInterface cachefreeLowerBound = new SimpleLowerBoundEstimationWithoutCaching();
		BermudanSwaption subBermudan = this.bermudanSwaption
				.getCloneWithModifiedStartingAndFinalPeriod(startingOptionPeriod,terminationOptionPeriod);
		subBermudan.setValuationMethod(cachefreeLowerBound);
		// calculate discount factor
		RandomVariable discountFactor = modelSubsimulation.getNumeraire(evaluationTime)
				.div(modelSubsimulation.getMonteCarloWeights(evaluationTime));
		// calculate discounted option value
		return subBermudan.getValue(evaluationTime, modelSubsimulation).div(discountFactor)
				.getAverage();
	}



	/**
	 * This method creates a new model with a changed starting and final date and
	 * initializes the LIBOR rates according to the given path of the input model.
	 * 
	 * @param model                LIBORModelMonteCarloSimulationModel to be used.
	 * @param startingPeriod       Period when the new model should start.
	 * @param endPeriod            Period when the new model should stop
	 *                             (inclusive).
	 * @param path                 Path of the original model to be used for
	 *                             starting values.
	 * @param pathsOfSubsimulation Number of paths for the new model.
	 * @return The model to be used for subsimulations.
	 * @throws CalculationException
	 */
	private LIBORModelMonteCarloSimulationModel createSubsimulationModelTerminating(
			LIBORModelMonteCarloSimulationModel model, int startingPeriod, int endPeriod, int path,
			int pathsOfSubsimulation) throws CalculationException {
		// shorten time discretization
		TimeDiscretization shortenedLiborDiscretization = new TimeDiscretizationFromArray(Arrays
				.copyOfRange(model.getLiborPeriodDiscretization().getAsDoubleArray(), startingPeriod, endPeriod + 1));
		// get initial LIBOR rates as forward curve starting point
		RandomVariable[] forwardCurveRandomVariable = model.getLIBORs(startingPeriod);
		ForwardCurve forwardCurve = createForwardCurveFromRealization (forwardCurveRandomVariable,shortenedLiborDiscretization, path);
		LIBORMarketModel liborMarketModel = liborModelCreator(shortenedLiborDiscretization,forwardCurve, model);		
		// adjust process if less factors needed
		int numberOfLibors = liborMarketModel.getNumberOfLibors();
		int numberOfFactors;
		if (numberOfLibors<model.getBrownianMotion().getNumberOfFactors())
			numberOfFactors = numberOfLibors;
		else
			numberOfFactors = model.getBrownianMotion().getNumberOfFactors();
		BrownianMotion brownianMotion = new net.finmath.montecarlo.BrownianMotionLazyInit(shortenedLiborDiscretization,
				numberOfFactors, pathsOfSubsimulation, 1234 + path);
		MonteCarloProcess process = new EulerSchemeFromProcessModel(brownianMotion,
				EulerSchemeFromProcessModel.Scheme.PREDICTOR_CORRECTOR);
		return new LIBORMonteCarloSimulationFromLIBORModel(liborMarketModel, process);
	}

	private LIBORMarketModel liborModelCreator(TimeDiscretization shortenedLiborDiscretization,
			ForwardCurve forwardCurve, LIBORModelMonteCarloSimulationModel model) throws CalculationException {
		DiscountCurveFromForwardCurve discountCurve = new DiscountCurveFromForwardCurve(forwardCurve);
		// create new model
		// Create a volatility structure v[i][j] = sigma_j(t_i)
		double[][] volatility = createVolatilityMatrix(0.3, 0.2, 0.25, forwardCurve, shortenedLiborDiscretization);
		LIBORVolatilityModelFromGivenMatrix volatilityModel = new LIBORVolatilityModelFromGivenMatrix(
				shortenedLiborDiscretization, shortenedLiborDiscretization, volatility);
		// Create a correlation model rho_{i,j} = exp(-a * abs(T_i-T_j))
		LIBORCorrelationModelExponentialDecay correlationModel = new LIBORCorrelationModelExponentialDecay(
				shortenedLiborDiscretization, shortenedLiborDiscretization, model.getNumberOfFactors(), 0.3);
		//  Combine volatility model and correlation model to a covariance model
		LIBORCovarianceModelFromVolatilityAndCorrelation covarianceModel = new LIBORCovarianceModelFromVolatilityAndCorrelation(
				shortenedLiborDiscretization, shortenedLiborDiscretization, volatilityModel, correlationModel);
		// create AAD random variable factory
		RandomVariableDifferentiableAADFactory randomVariableFactory = new RandomVariableDifferentiableAADFactory(
				new RandomVariableFactory());
		return new LIBORMarketModelFromCovarianceModel(shortenedLiborDiscretization, null,
				forwardCurve, discountCurve, randomVariableFactory, covarianceModel, null, null);
	}

	private double[][] createVolatilityMatrix(double a, double b, double c, ForwardCurve forwardCurve, TimeDiscretization shortenedLiborDiscretization) {
		int numberOfTimes = shortenedLiborDiscretization.getNumberOfTimes();
		int numberOfLIBORs = shortenedLiborDiscretization.getNumberOfTimes();
		double[][] volatility = new double[numberOfTimes][numberOfLIBORs];
		for (int timeIndex = 0; timeIndex < volatility.length; timeIndex++) {
			for (int liborIndex = 0; liborIndex < volatility[timeIndex].length; liborIndex++) {
				// Create a very simple volatility model here
				double time = shortenedLiborDiscretization.getTime(timeIndex);
				double maturity = shortenedLiborDiscretization.getTime(liborIndex);
				double timeToMaturity = maturity - time;
				double instVolatility;
				if (timeToMaturity <= 0) {
					instVolatility = 0; // This forward rate is already fixed, no volatility
				} else {
					instVolatility = (a + b * Math.exp(-c * timeToMaturity))
							* forwardCurve.getForward(null, shortenedLiborDiscretization.getTime(liborIndex)); // rescale
				}
				volatility[timeIndex][liborIndex] = instVolatility;
			}
		}
		return volatility;
	}
	private ForwardCurve createForwardCurveFromRealization(RandomVariable[] forwardCurveRandomVariable,
			TimeDiscretization shortenedLiborDiscretization, int path) {
		int numberOfShortendTimes = shortenedLiborDiscretization.getNumberOfTimes();
		double[] forwardArray = new double[numberOfShortendTimes];
		double periodLength = shortenedLiborDiscretization.getTimeStep(0);
		for (int i = 0; i < numberOfShortendTimes; i++)
			forwardArray[i] = forwardCurveRandomVariable[i].get(path);
		return ForwardCurveInterpolation.createForwardCurveFromForwards(
				"forwardCurve" /* name of the curve */, shortenedLiborDiscretization.getAsDoubleArray()
				/* fixings of the forward */, forwardArray, /* forwards */
				periodLength /* period lengths */
				);
	}
}
