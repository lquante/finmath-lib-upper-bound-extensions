package upperBoundMethods;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;

import bermudanSwaptionFramework.BermudanSwaption;
import bermudanSwaptionFramework.BermudanSwaptionValueEstimatorInterface;
import bermudanSwaptionFramework.SimplestExerciseStrategy;
import lowerBoundMethods.AbstractLowerBoundEstimationInputForUpperBound;
import lowerBoundMethods.AbstractLowerBoundEstimationWithoutCaching;
import lowerBoundMethods.SimpleLowerBoundEstimationWithoutCaching;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.montecarlo.BrownianMotion;
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
 * Implements the Andersen-Broadie method of martingale approximation using
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
	 * @param pathsSubsimulationsStepA number of subsimulation Paths in case 2a of
	 *                                 A-B algorithm, i.e. if exercise at current
	 *                                 simulation time.
	 * @param pathsSubsimulationsStepB number of subsimulation Paths in case 2b of
	 *                                 A-B algorithm, i.e. if no exercise at current
	 *                                 simulation time.
	 */

	public AndersenBroadieUpperBoundEstimation(AbstractLowerBoundEstimationInputForUpperBound lowerBoundMethod, double weightOfMartingale,
			int pathsSubsimulationsStepA, int pathsSubsimulationsStepB) {
		super(lowerBoundMethod,weightOfMartingale);
		this.pathsSubsimulationsStepA = pathsSubsimulationsStepA;
		this.pathsSubsimulationsStepB = pathsSubsimulationsStepB;

		// fixed exercise strategy, maybe add other implementations
		exerciseStrategy = new SimplestExerciseStrategy();
	}

	@Override
	protected double calculateMartingaleApproximation(int evaluationPeriod, LIBORModelMonteCarloSimulationModel model,
			RandomVariable[] cacheUnderlying, RandomVariable[] cacheOptionValues, RandomVariable[] triggerValues)
			throws CalculationException {
		int numberOfSimulations = cacheOptionValues[0].getRealizations().length;

		// determine number of martingale components to be estimated
		int numberOfOptionPeriods = this.bermudanSwaption.getFixingDates().length;

		// initialize martingale as lower bound value for period 0 and 1.
		ArrayList<RandomVariable> martingaleCache = new ArrayList<RandomVariable>();
		RandomVariable martingale = cacheOptionValues[evaluationPeriod];
		martingaleCache.add(martingale);
		if (evaluationPeriod + 1 < cacheOptionValues.length) {
			martingale = cacheOptionValues[evaluationPeriod + 1];
			martingaleCache.add(martingale);
		}

		for (int optionPeriod = 2 + evaluationPeriod; optionPeriod < numberOfOptionPeriods; optionPeriod++) {
			// initialize values for subsimulation
			// shifting if swaption starts in later period
			double currentFixingDate = this.bermudanSwaption.getFixingDates()[optionPeriod];
			int modelPeriod = model.getTimeIndex(currentFixingDate);

			int martingaleIndex = 2;

			// Arrays to store pathwise results to be transformed to RandomVariable
			double[] discountedExerciseValueArray = new double[numberOfSimulations];
			double[] discountedFutureExerciseValueArray = new double[numberOfSimulations];
			double[] previousExerciseIndicatorArray = new double[numberOfSimulations];

			// Parallelized implementation of pathwise part of A-B algorithm:
			int optionPeriodForInput = optionPeriod;
			int martingaleIndexForInput = martingaleIndex;
			IntStream.range(0, numberOfSimulations).parallel().forEach(path -> {
				double[] results = null;

				try {
					results = executeSubsimulation(model, path, optionPeriodForInput, modelPeriod, currentFixingDate,
							triggerValues, numberOfOptionPeriods, cacheUnderlying, martingaleIndexForInput);
				} catch (CalculationException e) {
					e.printStackTrace();
				}
				discountedExerciseValueArray[path] = results[0];
				discountedFutureExerciseValueArray[path] = results[1];
				previousExerciseIndicatorArray[path] = results[2];
			});

			RandomVariable discountedExerciseValue = new RandomVariableFromDoubleArray(currentFixingDate,
					discountedExerciseValueArray);
			RandomVariable discountedFutureExerciseValue = new RandomVariableFromDoubleArray(currentFixingDate,
					discountedFutureExerciseValueArray);
			RandomVariable previousExerciseIndicator = new RandomVariableFromDoubleArray(currentFixingDate,
					previousExerciseIndicatorArray);

			RandomVariable previousExerciseValue = cacheUnderlying[optionPeriod - 1];

			// calculate value of martingale extimation
			martingale = martingaleCache.get(martingaleIndex - 1).sub(cacheOptionValues[martingaleIndex - 1])
					.add(discountedExerciseValue)
					.sub((discountedFutureExerciseValue.add(previousExerciseValue)).mult(previousExerciseIndicator));
			martingaleCache.add(martingale);
			martingaleIndex += 1;

		}

		RandomVariable delta0 = model.getRandomVariableForConstant(0);
		// calculate the maximum of the estimated martingale for each remaining period
		for (int optionPeriod = 0 + evaluationPeriod; optionPeriod < numberOfOptionPeriods; optionPeriod++)
			delta0 = delta0.floor(cacheUnderlying[optionPeriod].sub(martingaleCache.get(optionPeriod)));

		// return average as the approximation of the martingale:
		return delta0.getAverage();

	}

	/**
	 * Method to be called by parallelized execution interface for subsimulations
	 * 
	 * @param model
	 * @param path
	 * @param optionPeriod
	 * @param modelPeriod
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
			int modelPeriod, double currentFixingDate, RandomVariable[] triggerValues, int numberOfOptionPeriods,
			RandomVariable[] cacheUnderlying, int martingaleIndex) throws CalculationException {
		// check trigger values for each path of the main simulation:

		double discountedExerciseValue = 0;
		double discountedFutureExerciseValue = 0;
		if (triggerValues[optionPeriod].get(path) >= 0) // case 2a
		{
			// create model for subsimulation
			LIBORModelMonteCarloSimulationModel modelStepA = createSubsimulationModel(model, modelPeriod, path,
					pathsSubsimulationsStepA);

			// create option with lower bound method without caching
			BermudanSwaptionValueEstimatorInterface cachefreeLowerBound = new SimpleLowerBoundEstimationWithoutCaching();
			BermudanSwaption bermudanA = this.bermudanSwaption.getCloneWithModifiedStartingPeriod(modelPeriod);
			bermudanA.setValuationMethod(cachefreeLowerBound);
			// calculate discount factor
			RandomVariable discountFactor = modelStepA.getNumeraire(currentFixingDate)
					.div(modelStepA.getMonteCarloWeights(currentFixingDate));
			// calculate discounted option value
			discountedExerciseValue = bermudanA.getValue(currentFixingDate, modelStepA).div(discountFactor)
					.getAverage();

		} else // case 2b
		{
			// create model terminating as soon as triggerValues[path+i] >=0 for
			// i>=optionPeriod+1
			// determine termination period:
			int terminationPeriod = optionPeriod;
			while (terminationPeriod < numberOfOptionPeriods - 1) {
				if (triggerValues[terminationPeriod].get(path) < 0)
					terminationPeriod += 1;
				else
					break;
			}
			double terminationTime = this.bermudanSwaption.getFixingDates()[terminationPeriod];
			int terminationTimeIndex = model.getTimeIndex(terminationTime);
			// calculate discounted option value:

			discountedExerciseValue = cacheUnderlying[optionPeriod].get(path);

			// create model for subsimulations (only if not only degenerated swaptions, thus
			// at least 2 additional periods necessary)
			if (optionPeriod + 2 < terminationPeriod && terminationPeriod < numberOfOptionPeriods) {
				LIBORModelMonteCarloSimulationModel modelStepB = createSubsimulationModelTerminating(model, modelPeriod,
						terminationTimeIndex, path, pathsSubsimulationsStepB);

				// create option
				double futureExerciseTime = modelStepB.getTime(1);

				BermudanSwaptionValueEstimatorInterface cachefreeLowerBound = new SimpleLowerBoundEstimationWithoutCaching();

				BermudanSwaption bermudanB = this.bermudanSwaption
						.getCloneWithModifiedStartingAndFinalPeriod(optionPeriod, terminationPeriod - 1);
				bermudanB.setValuationMethod(cachefreeLowerBound);
				// calculate discount factor
				RandomVariable discountFactor = modelStepB.getNumeraire(futureExerciseTime)
						.div(modelStepB.getMonteCarloWeights(futureExerciseTime));
				// calculate discounted option value

				discountedFutureExerciseValue = bermudanB.getValue(futureExerciseTime, modelStepB).div(discountFactor)
						.getAverage();
			}

		}
		// determine previous exercise indicator for martingale calculation

		double previousExerciseIndicator;
		if (triggerValues[martingaleIndex - 1].get(path) >= 0)
			previousExerciseIndicator = 1;
		else
			previousExerciseIndicator = 0;

		// return results as a triple
		double[] results = { discountedExerciseValue, discountedFutureExerciseValue, previousExerciseIndicator };
		return results;
	}

	/**
	 * This method creates a new model with a changed starting date and initializes
	 * the LIBOR rates according to the given path of the input model.
	 * 
	 * @param LIBORModelMonteCarloSimulationModel to be used
	 * @param startingPeriod                      Period when the new model should
	 *                                            start.
	 * @param path                                Path of the original model to be
	 *                                            used for starting values.
	 * @param pathsOfSubsimulation                Number of paths for the new model.
	 * @return The model to be used for subsimulations.
	 * @throws CalculationException
	 */
	private LIBORModelMonteCarloSimulationModel createSubsimulationModel(LIBORModelMonteCarloSimulationModel model,
			int startingPeriod, int path, int pathsOfSubsimulation) throws CalculationException {
		int numberOfTimes = model.getLiborPeriodDiscretization().getNumberOfTimes();

		return createSubsimulationModelTerminating(model, startingPeriod, numberOfTimes, path, pathsOfSubsimulation);
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
	public LIBORModelMonteCarloSimulationModel createSubsimulationModelTerminating(
			LIBORModelMonteCarloSimulationModel model, int startingPeriod, int endPeriod, int path,
			int pathsOfSubsimulation) throws CalculationException {

		TimeDiscretization shortenedLiborDiscretization = new TimeDiscretizationFromArray(Arrays
				.copyOfRange(model.getLiborPeriodDiscretization().getAsDoubleArray(), startingPeriod, endPeriod + 1));
		int numberOfShortendTimes = shortenedLiborDiscretization.getNumberOfTimes();
		// get initial LIBOR rates as forward curve starting point
		// get inital value of forward curve
		RandomVariable[] forwardCurveRandomVariable = model.getLIBORs(startingPeriod);
		double[] forwardArray = new double[numberOfShortendTimes];
		double periodLength = model.getLiborPeriodDiscretization().getTimeStep(startingPeriod);
		for (int i = 0; i < numberOfShortendTimes; i++)
			forwardArray[i] = forwardCurveRandomVariable[i].get(path);

		ForwardCurve forwardCurve = ForwardCurveInterpolation.createForwardCurveFromForwards(
				"forwardCurve" /* name of the curve */, shortenedLiborDiscretization.getAsDoubleArray()
				/* fixings of the forward */, forwardArray, /* forwards */
				periodLength /* period lengths */
		);

		DiscountCurveFromForwardCurve discountCurve = new DiscountCurveFromForwardCurve(forwardCurve);
		// create new model
		/*
		 * Create a volatility structure v[i][j] = sigma_j(t_i)
		 */
		double[][] volatility = new double[numberOfShortendTimes][numberOfShortendTimes];
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
					instVolatility = (0.3 + 0.2 * Math.exp(-0.25 * timeToMaturity))
							* forwardCurve.getForward(null, shortenedLiborDiscretization.getTime(liborIndex)); // rescale
					// values
				}

				// Store
				volatility[timeIndex][liborIndex] = instVolatility;
			}
		}
		LIBORVolatilityModelFromGivenMatrix volatilityModel = new LIBORVolatilityModelFromGivenMatrix(
				shortenedLiborDiscretization, shortenedLiborDiscretization, volatility);
		/*
		 * Create a correlation model rho_{i,j} = exp(-a * abs(T_i-T_j))
		 */
		LIBORCorrelationModelExponentialDecay correlationModel = new LIBORCorrelationModelExponentialDecay(
				shortenedLiborDiscretization, shortenedLiborDiscretization, model.getNumberOfFactors(), 0.3);

		/*
		 * Combine volatility model and correlation model to a covariance model
		 */
		LIBORCovarianceModelFromVolatilityAndCorrelation covarianceModel = new LIBORCovarianceModelFromVolatilityAndCorrelation(
				shortenedLiborDiscretization, shortenedLiborDiscretization, volatilityModel, correlationModel);

		// create AAD random variable factory
		RandomVariableDifferentiableAADFactory randomVariableFactory = new RandomVariableDifferentiableAADFactory(
				new RandomVariableFactory());

		LIBORMarketModel liborMarketModel = new LIBORMarketModelFromCovarianceModel(shortenedLiborDiscretization, null,
				forwardCurve, discountCurve, randomVariableFactory, covarianceModel, null, null);

		// adjust process
		int seed = 1234;
		BrownianMotion brownianMotion = new net.finmath.montecarlo.BrownianMotionLazyInit(shortenedLiborDiscretization,
				model.getBrownianMotion().getNumberOfFactors(), pathsOfSubsimulation, seed + path);
		MonteCarloProcess process = new EulerSchemeFromProcessModel(brownianMotion,
				EulerSchemeFromProcessModel.Scheme.PREDICTOR_CORRECTOR);

		return new LIBORMonteCarloSimulationFromLIBORModel(liborMarketModel, process);

	}
}
