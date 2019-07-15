package drafts;

import java.util.ArrayList;
import java.util.Arrays;


import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.montecarlo.BrownianMotion;
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

public class AndersenBroadieUpperBoundEstimation extends AbstractUpperBoundEstimation {
	AbstractSimpleBoundEstimation lowerBoundMethod;
	SimplestExerciseStrategy exerciseStrategy;
	private int numberOfSubsimulationsStepA;
	private int numberOfSubsimulationsStepB;
	private int numberOfSimulations;

	/**
	 * @param AbstractSimpleBoundEstimation lowerBoundMethod The lower bound method
	 *                                      to be used as a basis for the upper
	 *                                      bound.
	 */
	public AndersenBroadieUpperBoundEstimation(AbstractSimpleBoundEstimation lowerBoundMethod,
			int numberOfSubsimulationsStepA, int numberOfSubsimulationsStepB) {
		super(lowerBoundMethod);
		this.numberOfSubsimulationsStepA = numberOfSubsimulationsStepA;
		this.numberOfSubsimulationsStepB = numberOfSubsimulationsStepB;
		exerciseStrategy = new SimplestExerciseStrategy();
	}

	@Override
	protected RandomVariable calculateOptionValue(int period, LIBORModelMonteCarloSimulationModel model,
			RandomVariable[] cacheUnderlying, RandomVariable[] cacheOptionValues, RandomVariable[] triggerValues)
			throws CalculationException {
		numberOfSimulations = cacheOptionValues[0].getRealizations().length;
		
		double sumOfDeltas = 0;
		double delta0 = 0;
		for (int path = 0; path < numberOfSimulations; path++) {

			// initialize martingale as lower bound value for period 0 and 1.
			ArrayList<Double> martingaleCache = new ArrayList<Double>();
			double martingale = cacheOptionValues[period].get(path);
			martingaleCache.add(martingale);
			if (this.bermudanOption.getFixingDates().length > 1 && period+1<cacheOptionValues.length) {
				martingale = cacheOptionValues[period+1].get(path);
				martingaleCache.add(martingale);
			}

			// determine number of additional martingale components to be estimated
			int numberOfForwardPeriods = this.bermudanOption.getFixingDates().length - 1;

			// create model for subsimulations
			LIBORModelMonteCarloSimulationModel modelStepA = createSubsimulationModel(model, period, path,
					numberOfSubsimulationsStepA);
			
			
			
			for (int forwardPeriod = 2+period; forwardPeriod <= numberOfForwardPeriods; forwardPeriod++) {
				// initialize values for subsimulation
				double simulationTime = model.getTime(forwardPeriod);
				double discountedExerciseValue = 0;
				double discountedFutureExerciseValue = 0;
				int martingaleIndex = 2;
				// check trigger values for each path of the main simulation:
				if (triggerValues[period].get(path) >= 0) // case 2a
				{
					

					// create option
					BermudanSwaptionValueEstimatorInterface lowerBoundMethodA = new SimpleLowerBoundEstimation();
					BermudanSwaption BermudanA = this.bermudanOption
							.getBermudanSwaptionWithChangedValuationMethod(lowerBoundMethodA);
					BermudanSwaption restartedBermudanA = BermudanA.getCloneWithModifiedStartingPeriod(simulationTime);
					RandomVariable discountFactor = modelStepA.getNumeraire(simulationTime)
							.div(modelStepA.getMonteCarloWeights(forwardPeriod));
					// calculate discounted option value:

					discountedExerciseValue = restartedBermudanA.getValue(simulationTime, modelStepA).div(discountFactor)
							.getAverage();

				} else // case 2b
				{
					// create model terminating as soon as triggerValues[path+i] >=0 for i>=1
					int terminationPeriod;
					for (terminationPeriod = period;(triggerValues[terminationPeriod].get(path)<0);terminationPeriod++)
					{
						if(terminationPeriod>=numberOfForwardPeriods)
							break;
					}
					// calculate discounted option value:

					discountedExerciseValue = cacheUnderlying[forwardPeriod].get(path);

					
					// create model for subsimulations
					if (forwardPeriod<terminationPeriod)
					{
					LIBORModelMonteCarloSimulationModel modelStepB = createSubsimulationModelTerminating(model, period, terminationPeriod, path,
							numberOfSubsimulationsStepB);

						
						
					// create option
					BermudanSwaptionValueEstimatorInterface lowerBoundMethodB = new SimpleLowerBoundEstimation();
					BermudanSwaption BermudanB = this.bermudanOption
							.getBermudanSwaptionWithChangedValuationMethod(lowerBoundMethodB);
					BermudanSwaption restartedBermudanB = BermudanB.getCloneWithModifiedStartingPeriod(simulationTime);
					double forwardOptionTime = modelStepB.getTime(forwardPeriod-period);
					RandomVariable discountFactor = modelStepB.getRandomVariableForConstant(0);
					 discountFactor = modelStepB.getNumeraire(forwardOptionTime)
							.div(modelStepB.getMonteCarloWeights(forwardPeriod-period+1));
					
					
					// check if forwardPeriod<final period
					if (terminationPeriod < forwardPeriod & forwardPeriod+1 < numberOfForwardPeriods) {
						RandomVariable valueOptionB = restartedBermudanB.getValue(forwardOptionTime, modelStepB);
						discountedFutureExerciseValue = valueOptionB
								.div(discountFactor).getAverage();
					}
					}
				}

				int previousExerciseIndicator;
				if (triggerValues[martingaleIndex-1].get(path)>=0)
					previousExerciseIndicator = 1;
				else
					previousExerciseIndicator = 0;
				double previousExerciseValue = 0;
				martingale = martingaleCache.get(martingaleIndex-1) - (cacheOptionValues[martingaleIndex-1].get(path))
						+ (discountedExerciseValue)-discountedFutureExerciseValue-previousExerciseIndicator*previousExerciseValue;
				martingaleCache.add(martingale);
				martingaleIndex +=1;
				previousExerciseValue = discountedExerciseValue;
			}
			
			// calculate the maximum from the simulated martingale for each remaining period
			for (int forwardPeriod = 0; forwardPeriod <= numberOfForwardPeriods-period; forwardPeriod++)
				delta0 = Math.max(delta0,
						cacheUnderlying[forwardPeriod].get(path) - martingaleCache.get(forwardPeriod));
			sumOfDeltas += delta0;
			
		}

		double deltaZeroApproximation = sumOfDeltas / numberOfSimulations;
		// Note that values is a relative price - no numeraire division is required
		RandomVariable numeraireAtPeriod = model.getNumeraire(model.getTime(period));
		RandomVariable monteCarloProbabilitiesAtSimulationTime = model.getMonteCarloWeights(period);
		RandomVariable discountFactor = numeraireAtPeriod.div(monteCarloProbabilitiesAtSimulationTime);
		
		RandomVariable upperBoundValue = cacheOptionValues[period].mult(discountFactor).add(deltaZeroApproximation);

		return upperBoundValue;

	}

	private LIBORModelMonteCarloSimulationModel createSubsimulationModel(LIBORModelMonteCarloSimulationModel model,
			int period, int path, int numberOfSubsimulations) throws CalculationException {
				int numberOfTimes = model.getLiborPeriodDiscretization().getNumberOfTimes();
				
		return createSubsimulationModelTerminating( model,
			 period,numberOfTimes , path, numberOfSubsimulations);
	}

	public LIBORModelMonteCarloSimulationModel createSubsimulationModelTerminating(LIBORModelMonteCarloSimulationModel model,
			int startingPeriod, int endPeriod, int component, int numberOfSubsimulationPaths) throws CalculationException {

		
		TimeDiscretization shortenedLiborDiscretization = new TimeDiscretizationFromArray(Arrays
				.copyOfRange(model.getLiborPeriodDiscretization().getAsDoubleArray(), startingPeriod, endPeriod+1));
		int numberOfShortendTimes = shortenedLiborDiscretization.getNumberOfTimes();
		// get initial LIBOR rates as forward curve starting point
		// get inital value of forward curve
		RandomVariable[] forwardCurveRandomVariable = model.getLIBORs(startingPeriod);
		double[] forwardArray = new double[numberOfShortendTimes - 1];
		for (int i = 0; i < numberOfShortendTimes - 1; i++)
			forwardArray[i] = forwardCurveRandomVariable[i].get(component);

		ForwardCurve forwardCurve = ForwardCurveInterpolation.createForwardCurveFromForwards(
				"forwardCurve" /* name of the curve */,
				Arrays.copyOfRange(shortenedLiborDiscretization.getAsDoubleArray(), 1,
						numberOfShortendTimes)/* fixings of the forward */,
				forwardArray, /* forwards */
				model.getLiborPeriodDiscretization().getTimeStep(startingPeriod) /* period lengths */
		);

		DiscountCurveFromForwardCurve discountCurve = new DiscountCurveFromForwardCurve(forwardCurve);
		// create new model
		/*
		 * Create a volatility structure v[i][j] = sigma_j(t_i)
		 */
		double[][] volatility = new double[shortenedLiborDiscretization
				.getNumberOfTimeSteps()][shortenedLiborDiscretization.getNumberOfTimeSteps()];
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

		LIBORMarketModel liborMarketModel = new LIBORMarketModelFromCovarianceModel(shortenedLiborDiscretization, null,
				forwardCurve, discountCurve, covarianceModel, null, null);

		// adjust process
		BrownianMotion brownianMotion = new net.finmath.montecarlo.BrownianMotionLazyInit(
				shortenedLiborDiscretization, model.getBrownianMotion().getNumberOfFactors(),
				numberOfSubsimulationsStepA, 3141);
		MonteCarloProcess process = new EulerSchemeFromProcessModel(brownianMotion,
				EulerSchemeFromProcessModel.Scheme.PREDICTOR_CORRECTOR);

		return new LIBORMonteCarloSimulationFromLIBORModel(liborMarketModel, process);

	}
}
