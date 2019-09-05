package simulationMethods;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.Locale;

import bermudanswaptionframework.BermudanSwaption;
import bermudanswaptionframework.BermudanSwaptionValueEstimatorInterface;
import lowerboundmethods.SimpleLowerBoundEstimation;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import upperboundmethods.AndersenBroadieUpperBoundEstimation;
import upperboundmethods.DeltaHedgingUpperBound;

/**
 * Class to generate simulation factories.
 * @author Lennart Quante
 *
 */
public class SimulationFactory {
	// model parameters
	double lastTimePoint;
	double timeDiscretizationLength;
	double liborPeriodLength;
	// monte carlo parameters
	int numberOfPaths;
	int numberOfSubsimulationsStepA;
	int numberOfSubsimulationsStepB;
	// option parameters
	int numberOfExercisePeriods;
	double optionPeriodLength;
	
	
	// formatter for values

		private static DecimalFormat formatterTime = new DecimalFormat("00.00", new DecimalFormatSymbols(Locale.ENGLISH));
		private static DecimalFormat formatterValue = new DecimalFormat(" ##0.00000%;-##0.00000%",
				new DecimalFormatSymbols(Locale.ENGLISH));
		private static DecimalFormat formattterTime = new DecimalFormat(" ##0.000s;-##0.000",
				new DecimalFormatSymbols(Locale.ENGLISH));
		private static DecimalFormat formatterDeviation = new DecimalFormat(" 0.00000E00;-0.00000E00",
				new DecimalFormatSymbols(Locale.ENGLISH));

	/**
	 * Factory to provide simulation runners with the given parameters.
	 * 
	 * @param lastTimePoint
	 * @param timeDiscretizationLength
	 * @param liborPeriodLength
	 * @param numberOfPaths
	 * @param numberOfSubsimulationsStepA
	 * @param numberOfSubsimulationsStepB
	 * @param numberOfExercisePeriods
	 * @param optionPeriodLength
	 */
	public SimulationFactory(double lastTimePoint, double timeDiscretizationLength, double liborPeriodLength,
			int numberOfPaths, int numberOfSubsimulationsStepA, int numberOfSubsimulationsStepB,
			int numberOfExercisePeriods, double optionPeriodLength) {
		super();
		this.lastTimePoint = lastTimePoint;
		this.timeDiscretizationLength = timeDiscretizationLength;
		this.liborPeriodLength = liborPeriodLength;
		this.numberOfPaths = numberOfPaths;
		this.numberOfSubsimulationsStepA = numberOfSubsimulationsStepA;
		this.numberOfSubsimulationsStepB = numberOfSubsimulationsStepB;
		this.numberOfExercisePeriods = numberOfExercisePeriods;
		this.optionPeriodLength = optionPeriodLength;
	}

	public void testOptionFlatCurve(double rateShift, double swapRate) throws CalculationException {
		double[] forwardInterpolationRates = new double[] { swapRate - rateShift, swapRate - rateShift };
		double[] forwardInterpolationTimePoints = new double[] { 0, lastTimePoint };

		executeSimulationRunWithFlexibleRates(swapRate, forwardInterpolationRates, forwardInterpolationTimePoints);
	}

	private void executeSimulationRunWithFlexibleRates(double swaprate, double[] forwardInterpolationRates,
			double[] forwardInterpolationTimePoints) throws CalculationException {
		// set and print parameters
		TestModelFactory.setLastTime(lastTimePoint);
		TestModelFactory.setLiborRateTimeHorzion(lastTimePoint);
		TestModelFactory.setTimeDiscretizationPeriodLength(timeDiscretizationLength);
		TestModelFactory.setLiborPeriodLength(liborPeriodLength);
		TestModelFactory.setForwardInterpolationRates(forwardInterpolationRates);
		TestModelFactory.setForwardInterpolationTimePoints(forwardInterpolationTimePoints);
		TestModelFactory.setNumberOfPaths(numberOfPaths);
		System.out.println("Number of paths: " + numberOfPaths);
		System.out.println("Number of subsimulation paths step A: " + numberOfSubsimulationsStepA);
		System.out.println("Number of subsimulation paths step B: " + numberOfSubsimulationsStepB);
		System.out.println("Number of exercise periods: " + numberOfExercisePeriods);
		System.out.println("Time discretization period length: " + timeDiscretizationLength);
		System.out.println("LIBOR period length: " + liborPeriodLength);
		System.out.println("Initial forward rates: " + Arrays.toString(forwardInterpolationRates));
		System.out.println("Option period length: " + optionPeriodLength);
		System.out.println("Option swap rate: " + swaprate);
		executePrintSwaptionValuationMethods(swaprate); // run test
	}
	
	public void executePrintSwaptionValuationMethods(double swaprate) throws CalculationException {

		LIBORModelMonteCarloSimulationModel liborModel = TestModelFactory.createLIBORMarketModel();

		// set swaption parameters

		TestSwaptionFactory testSwaptionCreator = new TestSwaptionFactory(numberOfExercisePeriods, optionPeriodLength,
				swaprate);
		testSwaptionCreator.setNumberOfExercisePeriods(numberOfExercisePeriods);
		testSwaptionCreator.setPeriodLength(optionPeriodLength);
		testSwaptionCreator.setSwaprate(swaprate);

		// print head of comparison table
		System.out.println("Bermudan Swaption prices:\n");
		System.out.println(
				"FirstFixingDate\tLower Bound\tUpper Bound(AB)\t\tUpperBound(Deltas)\tDeviation(AB)\tDeviation(Delta subsimfree)");
		// "EvaluationDate Lower Bound Upper Bound(AB) Deviation(AB) ");

		for (int startIndexLIBOR = 1; startIndexLIBOR < liborModel.getNumberOfLibors()
				- numberOfExercisePeriods * optionPeriodLength / liborPeriodLength; startIndexLIBOR++) {
			double firstFixingDate = liborModel.getLiborPeriod(startIndexLIBOR);

			/*
			 * Value a bermudan swaption
			 */

			System.out.print(formatterTime.format(firstFixingDate) + "\t");

			SimpleLowerBoundEstimation lowerBound = new SimpleLowerBoundEstimation();
			BermudanSwaption testSwaption = testSwaptionCreator.constructBermudanSwaption(firstFixingDate, lowerBound);
			double lowerBoundValue = timingValuationTest(lowerBound, testSwaption, liborModel);
			// System.out.print(Arrays.toString(testSwaption.getExerciseProbabilities()));
			// AB upper bound approximation
			AndersenBroadieUpperBoundEstimation ABupperBound = new AndersenBroadieUpperBoundEstimation(lowerBound, 1,
					numberOfSubsimulationsStepA, numberOfSubsimulationsStepB);
			double ABupperBoundValue = timingValuationTest(ABupperBound, testSwaption, liborModel);
			// Upper bound using martingale construction via Delta hedging
			DeltaHedgingUpperBound DeltaUpperBound = new DeltaHedgingUpperBound(lowerBound, 1);
			double DeltaUpperBoundValue = timingValuationTest(DeltaUpperBound, testSwaption, liborModel);

			// Absolute duality gap AB method
			double deviationAB = Math.abs(lowerBoundValue - ABupperBoundValue);
			System.out.print(formatterDeviation.format(deviationAB) + "\t");

			// Absolute duality gap delta Hedge method
			double deviationDeltaHedge = Math.abs(lowerBoundValue - DeltaUpperBoundValue);
			System.out.println(formatterDeviation.format(deviationDeltaHedge) + "\t");

		}
	}

	// help method for calculation and printing of values and calculation time
	private static double timingValuationTest(BermudanSwaptionValueEstimatorInterface valuationMethod,
			BermudanSwaption swaption, LIBORModelMonteCarloSimulationModel liborModel) throws CalculationException {
		swaption.setValuationMethod(valuationMethod);
		long timingValuationStart = System.currentTimeMillis();
		double simulatedValue = swaption.getValue(liborModel);
		long timingValuationEnd = System.currentTimeMillis();
		double lastOperationTimingValuation = (timingValuationEnd - timingValuationStart) / 1000.0;
		System.out.print(formatterValue.format(simulatedValue) + "\t");
		System.out.print(formattterTime.format(lastOperationTimingValuation) + "\t");

		return simulatedValue;
	}

	


}
