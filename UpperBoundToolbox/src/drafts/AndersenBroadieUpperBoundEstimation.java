package drafts;

import java.util.ArrayList;
import java.util.Arrays;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.RandomVariableFactory;
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

public class AndersenBroadieUpperBoundEstimation extends AbstractUpperBoundEstimation {
	AbstractSimpleBoundEstimation lowerBoundMethod;
	SimplestExerciseStrategy exerciseStrategy;
	private int numberOfSubsimulationsStepA;
	private int numberOfSubsimulationsStepB;


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
	protected double calculateDeltaZero(int evaluationPeriod, LIBORModelMonteCarloSimulationModel model,
			RandomVariable[] cacheUnderlying, RandomVariable[] cacheOptionValues, RandomVariable[] triggerValues)
					throws CalculationException {
		int numberOfSimulations = cacheOptionValues[0].getRealizations().length;

		double sumOfDeltas = 0;
		double delta0 = 0;

		// pathwise calculation, as for every path depending on the exercise indicator step A or step B requires a different model (with the starting values of the path) and product
		// determine number of martingale components to be estimated
		int numberOfOptionPeriods = this.bermudanOption.getFixingDates().length;
		double startingPeriodBermudan = this.bermudanOption.getFixingDates()[0];
		int startingShift = model.getTimeIndex(startingPeriodBermudan);
		
		
		for (int path = 0; path < numberOfSimulations; path++) {

			// initialize martingale as lower bound value for period 0 and 1.
			ArrayList<Double> martingaleCache = new ArrayList<Double>();
			double martingale = cacheOptionValues[evaluationPeriod].get(path);
			martingaleCache.add(martingale);
			if (evaluationPeriod + 1 < cacheOptionValues.length) {
				martingale = cacheOptionValues[evaluationPeriod + 1].get(path);
				martingaleCache.add(martingale);
			}

		

			for (int optionPeriod = 2 + evaluationPeriod; optionPeriod < numberOfOptionPeriods; optionPeriod++) {
				// initialize values for subsimulation
				int modelPeriod = optionPeriod + startingShift; // shifting if swaption starts in later period
				double currentFixingDate = this.bermudanOption.getFixingDates()[optionPeriod];
				double discountedExerciseValue = 0;
				double discountedFutureExerciseValue = 0;
				int martingaleIndex = 2;
				// check trigger values for each path of the main simulation:
				if (triggerValues[optionPeriod].get(path) >= 0) // case 2a
				{
					// create model for subsimulation
					LIBORModelMonteCarloSimulationModel modelStepA = createSubsimulationModel(model,
							modelPeriod, path, numberOfSubsimulationsStepA);

					// create option
					BermudanSwaption bermudanA = this.bermudanOption.getCloneWithModifiedStartingPeriod(modelPeriod);
					// calculate discount factor
					RandomVariable discountFactor = modelStepA.getNumeraire(currentFixingDate)
							.div(modelStepA.getMonteCarloWeights(1));
					// calculate discounted option value
					discountedExerciseValue = bermudanA.getValue(currentFixingDate, modelStepA).div(discountFactor)
							.getAverage();

				} else // case 2b
				{
					// create model terminating as soon as triggerValues[path+i] >=0 for i>=1
					// determine termination period:
					int terminationPeriod= optionPeriod+1;
					while(terminationPeriod < numberOfOptionPeriods-1)
					{
						if(triggerValues[terminationPeriod].get(path) < 0)
							terminationPeriod+=1;
						else
							break;
					}
					
					// calculate discounted option value:

					discountedExerciseValue = cacheUnderlying[optionPeriod].get(path);

					// create model for subsimulations (only if not only degenerated swaptions, thus at least 2 additional periods necessary)
					if (optionPeriod+2< terminationPeriod && terminationPeriod < numberOfOptionPeriods) {
						LIBORModelMonteCarloSimulationModel modelStepB = createSubsimulationModelTerminating(model,
								modelPeriod, terminationPeriod+startingShift, path, numberOfSubsimulationsStepB);

						// create option
						double futureExerciseTime = modelStepB.getTime(1);
						BermudanSwaption bermudanB = this.bermudanOption
							.getCloneWithModifiedStartingAndFinalPeriod(optionPeriod, terminationPeriod-1);
						// calculate discount factor
						RandomVariable discountFactor = modelStepB.getNumeraire(futureExerciseTime)
								.div(modelStepB.getMonteCarloWeights(futureExerciseTime));
						// calculate discounted option value
						discountedFutureExerciseValue = bermudanB.getValue(futureExerciseTime, modelStepB).div(discountFactor).getAverage();						}
				}


				// determine previous exercise value for martingale calculation
				int previousExerciseIndicator;
				if (triggerValues[martingaleIndex - 1].get(path) >= 0)
					previousExerciseIndicator = 1;
				else
					previousExerciseIndicator = 0;
				double previousExerciseValue = cacheUnderlying[optionPeriod - 1].get(path);
				
				// calculate value of martingale extimation
				martingale = martingaleCache.get(martingaleIndex - 1)
						- (cacheOptionValues[martingaleIndex - 1].get(path)) + discountedExerciseValue
						-  previousExerciseIndicator * (discountedFutureExerciseValue+previousExerciseValue);
				martingaleCache.add(martingale);
				martingaleIndex += 1;

			}

			// calculate the maximum using the estimated martingale for each remaining period
			for (int forwardPeriod = 0; forwardPeriod < numberOfOptionPeriods - evaluationPeriod; forwardPeriod++)
				delta0 = Math.max(delta0,
						cacheUnderlying[forwardPeriod].get(path) - martingaleCache.get(forwardPeriod));
			sumOfDeltas += delta0;

		}
		// return average of the approximation of martingale:
		return sumOfDeltas / numberOfSimulations;

		

	}

	private LIBORModelMonteCarloSimulationModel createSubsimulationModel(LIBORModelMonteCarloSimulationModel model,
			int period, int path, int numberOfSubsimulations) throws CalculationException {
		int numberOfTimes = model.getLiborPeriodDiscretization().getNumberOfTimes();

		return createSubsimulationModelTerminating(model, period, numberOfTimes, path, numberOfSubsimulations);
	}

	public LIBORModelMonteCarloSimulationModel createSubsimulationModelTerminating(
			LIBORModelMonteCarloSimulationModel model, int startingPeriod, int endPeriod, int component,
			int numberOfSubsimulationPaths) throws CalculationException {

		TimeDiscretization shortenedLiborDiscretization = new TimeDiscretizationFromArray(Arrays
				.copyOfRange(model.getLiborPeriodDiscretization().getAsDoubleArray(), startingPeriod, endPeriod + 1));
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

		// create aad random variable factory
		RandomVariableDifferentiableAADFactory randomVariableFactory = new RandomVariableDifferentiableAADFactory(
				new RandomVariableFactory());

		LIBORMarketModel liborMarketModel = new LIBORMarketModelFromCovarianceModel(shortenedLiborDiscretization, null,
				forwardCurve, discountCurve, randomVariableFactory, covarianceModel, null, null);

		// adjust process
		int seed=1234;
		BrownianMotion brownianMotion = new net.finmath.montecarlo.BrownianMotionLazyInit(shortenedLiborDiscretization,
				model.getBrownianMotion().getNumberOfFactors(), numberOfSubsimulationPaths, seed+component);
		MonteCarloProcess process = new EulerSchemeFromProcessModel(brownianMotion,
				EulerSchemeFromProcessModel.Scheme.PREDICTOR_CORRECTOR);

		return new LIBORMonteCarloSimulationFromLIBORModel(liborMarketModel, process);

	}
}
