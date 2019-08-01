package tests;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAADFactory;
import net.finmath.montecarlo.interestrate.CalibrationProduct;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.LIBORMonteCarloSimulationFromLIBORModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelFromGivenMatrix;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel.Scheme;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

public class CreateTestModel {

	// all model parameters are set as fields and used in the creation method to
	// facilitate adjustments

	LIBORModelMonteCarloSimulationModel liborModel;
	// monte carlo simulation parameters
	static int numberOfPaths = 1000;
	static int numberOfFactors = 3; // PCA number of factors
	static int seed = 3141; // seed for stochastic driver

	static AbstractRandomVariableFactory randomVariableFactory = new RandomVariableDifferentiableAADFactory(
			new RandomVariableFactory());
	// state space
	static String stateSpace = LIBORMarketModelFromCovarianceModel.StateSpace.NORMAL.name();
	// simulation scheme
	private static Scheme scheme = EulerSchemeFromProcessModel.Scheme.PREDICTOR_CORRECTOR;

	// time tiscretization properties

	static double lastTime = 20.0;
	static double dt = 0.5;

	static // properties of the libor model
	double liborRateTimeHorzion = 20.0;
	static double liborPeriodLength = 0.5;
	final static double correlationDecayParam = 0.3;
	// points for interpolaton forward curve
	static double[] forwardInterpolationTimePoints = new double[] { 0.5, 1.0, 2.0, 5.0, 40.0 };
	static double[] forwardInterpolationRates = new double[] { 0.02, 0.02, 0.02, 0.02, 0.02 };
	// parameters of volatility model (instVolatility = (a + b * Math.exp(-c *
	// timeToMaturity)
	static double a = 0.3;
	static double b = 0.2;
	static double c = 0.25;
	// measure
	static String measure = LIBORMarketModelFromCovarianceModel.Measure.SPOT.name();

	public CreateTestModel() throws CalculationException {
		liborModel = createLIBORMarketModel();
	}

	@Test
	public void testModel() {

		// print some model properties
		TimeDiscretization timeDiscretization = liborModel.getTimeDiscretization();
		TimeDiscretization liborDiscretization = liborModel.getLiborPeriodDiscretization();
		System.out.println("TimeDiscretization:" + Arrays.toString(timeDiscretization.getAsDoubleArray()));
		System.out.println("LiborDiscretization:" + Arrays.toString(liborDiscretization.getAsDoubleArray()));

		// try to calculate all LIBOR rates
		for (int i = 0; i < liborDiscretization.getNumberOfTimeSteps(); i++)
			try {
				RandomVariable[] libors = liborModel.getLIBORs(i);
				System.out.println(i + "thLIBOR spot rate " + libors[i].getAverage());
			} catch (CalculationException e) {
				System.out.println("calculation of libors failed at timepoint" + liborDiscretization.getTime(i));
				e.printStackTrace();
			}
	}

	public static LIBORModelMonteCarloSimulationModel createLIBORMarketModel() throws CalculationException {

		/*
		 * Create the libor tenor structure and the initial values
		 */

		TimeDiscretizationFromArray liborPeriodDiscretization = new TimeDiscretizationFromArray(0.0,
				(int) (liborRateTimeHorzion / liborPeriodLength), liborPeriodLength);

		// Create the forward curve (initial value of the LIBOR market model)
		ForwardCurveInterpolation forwardCurveInterpolation = ForwardCurveInterpolation.createForwardCurveFromForwards(
				"forwardCurve" /* name of the curve */, forwardInterpolationTimePoints, forwardInterpolationRates,
				liborPeriodLength);

		/*
		 * Create a simulation time discretization
		 */

		TimeDiscretizationFromArray timeDiscretizationFromArray = new TimeDiscretizationFromArray(0.0,
				(int) (lastTime / dt), dt);

		/*
		 * Create a volatility structure v[i][j] = sigma_j(t_i)
		 */
		double[][] volatility = new double[timeDiscretizationFromArray.getNumberOfTimeSteps()][liborPeriodDiscretization
				.getNumberOfTimeSteps()];
		for (int timeIndex = 0; timeIndex < volatility.length; timeIndex++) {
			for (int liborIndex = 0; liborIndex < volatility[timeIndex].length; liborIndex++) {
				// Create a very simple volatility model here
				double time = timeDiscretizationFromArray.getTime(timeIndex);
				double maturity = liborPeriodDiscretization.getTime(liborIndex);
				double timeToMaturity = maturity - time;

				double instVolatility;
				if (timeToMaturity <= 0) {
					instVolatility = 0; // This forward rate is already fixed, no volatility
				} else {
					instVolatility = (a + b * Math.exp(-c * timeToMaturity))
							* forwardCurveInterpolation.getForward(null, liborPeriodDiscretization.getTime(liborIndex)); // rescale
					// values
				}

				// Store
				volatility[timeIndex][liborIndex] = instVolatility;
			}
		}
		LIBORVolatilityModelFromGivenMatrix volatilityModel = new LIBORVolatilityModelFromGivenMatrix(
				timeDiscretizationFromArray, liborPeriodDiscretization, volatility);
		/*
		 * Create a correlation model rho_{i,j} = exp(-a * abs(T_i-T_j))
		 */
		LIBORCorrelationModelExponentialDecay correlationModel = new LIBORCorrelationModelExponentialDecay(
				timeDiscretizationFromArray, liborPeriodDiscretization, numberOfFactors, correlationDecayParam);

		/*
		 * Combine volatility model and correlation model to a covariance model
		 */
		LIBORCovarianceModelFromVolatilityAndCorrelation covarianceModel = new LIBORCovarianceModelFromVolatilityAndCorrelation(
				timeDiscretizationFromArray, liborPeriodDiscretization, volatilityModel, correlationModel);

		// Set model properties
		Map<String, String> properties = new HashMap<>();

		// Choose the simulation measure
		properties.put("measure", measure);

		// Choose log normal model
		properties.put("stateSpace", stateSpace);

		// Empty array of calibration items - hence, model will use given covariance
		CalibrationProduct[] calibrationItems = new CalibrationProduct[0];

		/*
		 * Create corresponding LIBOR Market Model
		 */
		LIBORMarketModel liborMarketModel = new LIBORMarketModelFromCovarianceModel(liborPeriodDiscretization, null,
				forwardCurveInterpolation, new DiscountCurveFromForwardCurve(forwardCurveInterpolation),
				randomVariableFactory, covarianceModel, calibrationItems, properties);

		BrownianMotion brownianMotion = new net.finmath.montecarlo.BrownianMotionLazyInit(timeDiscretizationFromArray,
				numberOfFactors, numberOfPaths, seed /* seed */);

		EulerSchemeFromProcessModel process = new EulerSchemeFromProcessModel(brownianMotion, scheme);

		return new LIBORMonteCarloSimulationFromLIBORModel(liborMarketModel, process);
	}

	/**
	 * @return the numberOfPaths
	 */
	public static int getNumberOfPaths() {
		return numberOfPaths;
	}

	/**
	 * @param numberOfPaths the numberOfPaths to set
	 */
	public static void setNumberOfPaths(int numberOfPaths) {
		CreateTestModel.numberOfPaths = numberOfPaths;
	}

	/**
	 * @return the numberOfFactors
	 */
	public static int getNumberOfFactors() {
		return numberOfFactors;
	}

	/**
	 * @param numberOfFactors the numberOfFactors to set
	 */
	public static void setNumberOfFactors(int numberOfFactors) {
		CreateTestModel.numberOfFactors = numberOfFactors;
	}

	/**
	 * @return the seed
	 */
	public static int getSeed() {
		return seed;
	}

	/**
	 * @param seed the seed to set
	 */
	public static void setSeed(int seed) {
		CreateTestModel.seed = seed;
	}
}
