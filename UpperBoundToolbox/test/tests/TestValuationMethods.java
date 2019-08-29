package tests;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.Locale;

import org.junit.Assert;
import org.junit.Test;

import bermudanSwaptionFramework.BermudanSwaption;
import bermudanSwaptionFramework.BermudanSwaptionValueEstimatorInterface;
import lowerBoundMethods.SimpleLowerBoundEstimation;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import upperBoundMethods.AndersenBroadieUpperBoundEstimation;
import upperBoundMethods.DeltaHedgingUpperBound;

public class TestValuationMethods {
	// formatter for values

	private static DecimalFormat formatterTime = new DecimalFormat("00.00", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterValue = new DecimalFormat(" ##0.00000%;-##0.00000%",
			new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formattterTime = new DecimalFormat(" ##0.000s;-##0.000",
			new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterDeviation = new DecimalFormat(" 0.00000E00;-0.00000E00",
			new DecimalFormatSymbols(Locale.ENGLISH));

	// set tolerance for difference between upper and lower bound methods
	static double tolerance = 1; // should be tightened pending further improvement

	private static int numberOfPaths = 1000;
	private static double timeDiscretizationLength=0.25;
	private static double liborPeriodLength = 0.5;
	
	static int numberOfExercisePeriods = 4;
	private static int numberOfSubsimulationsStepA = 100;
	private static int numberOfSubsimulationsStepB = 100;
	private static double optionPeriodLength=1;

	

	@Test
	public void testSwaptionValuationMethods() throws CalculationException {
		// set parameters
		CreateTestModel.setNumberOfPaths(numberOfPaths);
		CreateTestModel.setTimeDiscretizationPeriodLength(timeDiscretizationLength);
		CreateTestModel.setLiborPeriodLength(liborPeriodLength);
		CreateTestBermudanSwaption.setNumberOfExercisePeriods(numberOfExercisePeriods);
		CreateTestBermudanSwaption.setPeriodLength(optionPeriodLength);
		TestValuationMethods.setNumberOfSubsimulationsStepA(numberOfSubsimulationsStepA);
		TestValuationMethods.setNumberOfSubsimulationsStepB(numberOfSubsimulationsStepB);
		
		
		executePrintSwaptionValuationMethods();
	}

	public static void executePrintSwaptionValuationMethods() throws CalculationException {

		LIBORModelMonteCarloSimulationModel liborModel = CreateTestModel.createLIBORMarketModel();
		
		// print head of comparison table
		System.out.println("Bermudan Swaption prices:\n");
		System.out.println(
				"FirstFixingDate\tLower Bound\tUpper Bound(AB)\t\tUpperBound(Deltas)\tDeviation(AB)\tDeviation(Delta subsimfree)");
		// "EvaluationDate Lower Bound Upper Bound(AB) Deviation(AB) ");

		for (int modelTimeIndexToStartOption = 1; modelTimeIndexToStartOption < liborModel.getNumberOfLibors()
				-numberOfExercisePeriods*optionPeriodLength/liborPeriodLength ; modelTimeIndexToStartOption++) {
			double firstFixingDate = liborModel.getLiborPeriod(modelTimeIndexToStartOption);
			CreateTestBermudanSwaption.setFirstFixingDate(firstFixingDate);
			CreateTestBermudanSwaption.setNumberOfExercisePeriods(numberOfExercisePeriods);
			BermudanSwaption testSwaption = CreateTestBermudanSwaption.createBermudanSwaption();

			/*
			 * Value a bermudan swaption
			 */

			System.out.print(formatterTime.format(firstFixingDate) + "\t");

			SimpleLowerBoundEstimation lowerBound = new SimpleLowerBoundEstimation();

			double lowerBoundValue = timingValuationTest(lowerBound, testSwaption, liborModel);
			System.out.print(Arrays.toString(testSwaption.getExerciseProbabilities()));
			// AB upper bound approximation
			AndersenBroadieUpperBoundEstimation ABupperBound = new AndersenBroadieUpperBoundEstimation(lowerBound,
					1,numberOfSubsimulationsStepA, numberOfSubsimulationsStepB);
			double ABupperBoundValue = timingValuationTest(ABupperBound, testSwaption, liborModel);
			// Upper bound using martingale construction via Delta hedging
			DeltaHedgingUpperBound DeltaUpperBound = new DeltaHedgingUpperBound(lowerBound,1);
			double DeltaUpperBoundValue = timingValuationTest(DeltaUpperBound, testSwaption, liborModel);

			// Absolute error AB method
			double deviationAB = Math.abs(lowerBoundValue - ABupperBoundValue);
			System.out.print(formatterDeviation.format(deviationAB) + "\t");

			// Absolute error delta Hedge method
			double deviationDeltaHedge = Math.abs(lowerBoundValue - DeltaUpperBoundValue);
			System.out.println(formatterDeviation.format(deviationDeltaHedge) + "\t");

			// Assert deviations
			Assert.assertEquals(lowerBoundValue, ABupperBoundValue, tolerance);
			Assert.assertEquals(lowerBoundValue, DeltaUpperBoundValue, tolerance);
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

	/**
	 * @return the numberOfSubsimulationsStepA
	 */
	public static int getNumberOfSubsimulationsStepA() {
		return numberOfSubsimulationsStepA;
	}

	/**
	 * @param numberOfSubsimulationsStepA the numberOfSubsimulationsStepA to set
	 */
	public static void setNumberOfSubsimulationsStepA(int numberOfSubsimulationsStepA) {
		TestValuationMethods.numberOfSubsimulationsStepA = numberOfSubsimulationsStepA;
	}

	/**
	 * @return the numberOfSubsimulationsStepB
	 */
	public static int getNumberOfSubsimulationsStepB() {
		return numberOfSubsimulationsStepB;
	}

	/**
	 * @param numberOfSubsimulationsStepB the numberOfSubsimulationsStepB to set
	 */
	public static void setNumberOfSubsimulationsStepB(int numberOfSubsimulationsStepB) {
		TestValuationMethods.numberOfSubsimulationsStepB = numberOfSubsimulationsStepB;
	}

}
