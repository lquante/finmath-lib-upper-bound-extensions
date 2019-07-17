/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 28 Feb 2015
 */
package tests;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

import org.junit.Assert;
import org.junit.Test;

import drafts.AndersenBroadieUpperBoundEstimation;
import drafts.BermudanSwaption;
import drafts.SimpleLowerBoundEstimation;
import modifiedFinmathFilesProcess.LinearAlgebra;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.interestrate.CalibrationProduct;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.LIBORMonteCarloSimulationFromLIBORModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelFromGivenMatrix;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.time.TimeDiscretizationFromArray;

/**
 * This class tests the LIBOR market model and products.
 *
 * @author Christian Fries
 * @author Lorenzo Torricelli
 */
public class ValuationOfBermudanSwaptionRudimentalTest {

	private static DecimalFormat formatterMaturity = new DecimalFormat("00.00",
			new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterValue = new DecimalFormat(" ##0.00000%;-##0.00000%",
			new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterDeviation = new DecimalFormat(" 0.00000E00;-0.00000E00",
			new DecimalFormatSymbols(Locale.ENGLISH));

	private final LIBORModelMonteCarloSimulationModel liborModel;

	private final int numberOfPaths = 1000;
	private final int numberOfSubsimulationsStepA = 100;
	private final int numberOfSubsimulationsStepB = 100;
	private final int numberOfFactors = 3; // PCA number of factors
	private final double correlationDecayParam = 0.3;
	private AbstractRandomVariableFactory randomVariableFactory= new RandomVariableFactory();

	public ValuationOfBermudanSwaptionRudimentalTest() throws CalculationException {
		liborModel = ValuationOfBermudanSwaptionRudimentalTest.createLIBORMarketModel(randomVariableFactory, numberOfPaths, numberOfFactors,
				correlationDecayParam);
		
	}

	@Test
	public void testSwaption() throws CalculationException {
		// set tolerance for difference upper and lower bound
		double tolerance = 1; // should be tightened pending further improvement
		
		//print head of comparison table
		System.out.println("Bermudan Swaption prices:\n");
		System.out.println(
				//"EvaluationDate      Lower Bound       Upper Bound(AB)        Upper Bound(subsimfree)          Deviation(AB)        Deviation(subsimfree)");
				"EvaluationDate      Lower Bound       Upper Bound(AB)                  Deviation(AB)        ");
		int numberOfPeriods = 10;
		// Create libor Market model
		for (int maturityIndex = 1; maturityIndex < liborModel.getNumberOfLibors() - numberOfPeriods; maturityIndex++) {
			double exerciseDate = liborModel.getLiborPeriod(maturityIndex);
			
			
			

			// Create a rudimental bermudan swaption

			double[] fixingDates = new double[numberOfPeriods];
			double[] paymentDates = new double[numberOfPeriods]; // simply a merge of fixing and payment dates (which
																	// obviously overlap)

			double[] notionals = new double[numberOfPeriods + 1];
			boolean[] startingDate = new boolean[numberOfPeriods + 1];
			double[] swapTenor = new double[numberOfPeriods]; // to be passed to the analytical approximation method
			double swapPeriodLength = 0.5; // Use instead -> liborMarketModel.getLiborTimediscretisation to make the
											// libor discretisation
			// coincide with the swap time discretisation
			
			
			
			for (int periodStartIndex = 0; periodStartIndex < numberOfPeriods; periodStartIndex++) {
				fixingDates[periodStartIndex] = exerciseDate + periodStartIndex * swapPeriodLength;
				paymentDates[periodStartIndex] = exerciseDate + (periodStartIndex + 1) * swapPeriodLength;
				swapTenor[periodStartIndex] = swapPeriodLength;
				notionals[periodStartIndex] = 1;
				startingDate[periodStartIndex] = true;
			}
			notionals[numberOfPeriods] = 1;
			startingDate[numberOfPeriods] = true;
			// swapTenor[numberOfPeriods] = exerciseDate + numberOfPeriods *
			// swapPeriodLength;

			// Swaptions swap rate
			double swaprate = 0.02; // getParSwaprate(liborModel, swapTenor);

			// Set swap rates for each period
			double[] swaprates = new double[numberOfPeriods];
			for (int periodStartIndex = 0; periodStartIndex < numberOfPeriods; periodStartIndex++) {
				swaprates[periodStartIndex] = swaprate;
			}

			/*
			 * Value a bermudan swaption
			 */
			
			System.out.print(formatterMaturity.format(exerciseDate) + "          ");
			// Value with lower bound method

			// change regression basis functions via enum:

			// choose actual basis function method here: BasisFunctionType.
			// default is using forward rates

			// SwapRates
			// ForwardRates (with cross)

			// more types to be developed
		
			SimpleLowerBoundEstimation lowerBound = new SimpleLowerBoundEstimation();

			String currency = "EURO";
			BermudanSwaption swaptionMonteCarlo = new BermudanSwaption(currency, startingDate, fixingDates,
					swapTenor, paymentDates, notionals, true, swaprates, lowerBound);
			double lowerBoundValue = swaptionMonteCarlo.getValue(liborModel);
			System.out.print(formatterValue.format(lowerBoundValue) + "          ");
			
			// print exercise probabilities
			//System.out.println("Exercise probabilites:"+Arrays.toString(lowerBound.getExerciseProbablities()));
			
			SimpleLowerBoundEstimation lowerBoundForUpper = new SimpleLowerBoundEstimation();
			  // Value of the AB upper bound approximation
			BermudanSwaption upperBoundSwaptionMonteCarlo = new BermudanSwaption(currency, startingDate, fixingDates,
					swapTenor, paymentDates, notionals, true, swaprates, lowerBoundForUpper);
			AndersenBroadieUpperBoundEstimation ABupperBound = new AndersenBroadieUpperBoundEstimation(
					lowerBoundForUpper, numberOfSubsimulationsStepA, numberOfSubsimulationsStepB);
			upperBoundSwaptionMonteCarlo.setValuationMethod(ABupperBound);
			double ABupperBoundValue = upperBoundSwaptionMonteCarlo.getValue(liborModel);
			System.out.print(formatterValue.format(ABupperBoundValue) + "          ");

//			// Value of the sub simulation free upper bound approximation
//			SubSimulationFreeUpperBoundEstimation subSimFreeUpperBound = new SubSimulationFreeUpperBoundEstimation(lowerBound, 100, 100);
//			swaptionMonteCarlo.setValuationMethod(subSimFreeUpperBound);
//			double subSimFreeUpperBoundValue = swaptionMonteCarlo.getValue(liborModel);
//			System.out.print(formatterValue.format(subSimFreeUpperBoundValue) + "          ");
//			
	
			// Absolute error AB method
			double deviationAB = Math.abs(lowerBoundValue - ABupperBoundValue);
			System.out.println(formatterDeviation.format(deviationAB) + "          ");

//			// Absolute error subsimfree method
//			double deviationSubsimfree = Math.abs(lowerBoundValue - subSimFreeUpperBoundValue);
//			System.out.println(formatterDeviation.format(deviationSubsimfree) + "          ");
//			
			Assert.assertEquals(lowerBoundValue, ABupperBoundValue, tolerance);
		}
	}

	public static LIBORModelMonteCarloSimulationModel createLIBORMarketModel(AbstractRandomVariableFactory randomVariableFactory,int numberOfPaths, int numberOfFactors,
			double correlationDecayParam) throws CalculationException {

		/*
		 * Create the libor tenor structure and the initial values
		 */
		double liborPeriodLength = 0.5;
		double liborRateTimeHorzion = 20.0;
		TimeDiscretizationFromArray liborPeriodDiscretization = new TimeDiscretizationFromArray(0.0,
				(int) (liborRateTimeHorzion / liborPeriodLength), liborPeriodLength);

		// Create the forward curve (initial value of the LIBOR market model)
		ForwardCurveInterpolation forwardCurveInterpolation = ForwardCurveInterpolation.createForwardCurveFromForwards(
				"forwardCurve" /* name of the curve */,
				new double[] { 0.5, 1.0, 2.0, 5.0, 40.0 } /* fixings of the forward */,
				new double[] { 0.02, 0.02, 0.02, 0.02, 0.02 } /* forwards */,
				liborPeriodLength /* tenor / period length */
		);

		/*
		 * Create a simulation time discretization
		 */
		double lastTime = 20.0;
		double dt = 0.5;

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
					instVolatility = (0.3 + 0.2 * Math.exp(-0.25 * timeToMaturity))
							* forwardCurveInterpolation.getForward(null, liborPeriodDiscretization.getTime(liborIndex)); // rescale
																								// values
				}

				// Store
				volatility[timeIndex][liborIndex] = instVolatility;
			}
		}
		LIBORVolatilityModelFromGivenMatrix volatilityModel = new LIBORVolatilityModelFromGivenMatrix(
				timeDiscretizationFromArray, liborPeriodDiscretization, volatility);
		// LIBORVolatilityModel volatilityModel = new
		// LIBORVolatilityModelFourParameterExponentialForm(timeDiscretizationFromArray,
		// liborPeriodDiscretization, 0.00, 0.0, 0.0, 0.01, false);

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
		properties.put("measure", LIBORMarketModelFromCovarianceModel.Measure.SPOT.name());

		// Choose log normal model
		properties.put("stateSpace", LIBORMarketModelFromCovarianceModel.StateSpace.NORMAL.name());

		// Empty array of calibration items - hence, model will use given covariance
		CalibrationProduct[] calibrationItems = new CalibrationProduct[0];

		/*
		 * Create corresponding LIBOR Market Model
		 */
		LIBORMarketModel liborMarketModel = new LIBORMarketModelFromCovarianceModel(liborPeriodDiscretization,
				null, forwardCurveInterpolation, new DiscountCurveFromForwardCurve(forwardCurveInterpolation),randomVariableFactory,
				covarianceModel, calibrationItems, properties);

		BrownianMotion brownianMotion = new net.finmath.montecarlo.BrownianMotionLazyInit(timeDiscretizationFromArray,
				numberOfFactors, numberOfPaths, 3141 /* seed */);

		EulerSchemeFromProcessModel process = new EulerSchemeFromProcessModel(brownianMotion,
				EulerSchemeFromProcessModel.Scheme.PREDICTOR_CORRECTOR);

		return new LIBORMonteCarloSimulationFromLIBORModel(liborMarketModel, process);
	}
}
