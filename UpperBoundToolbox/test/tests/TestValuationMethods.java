package tests;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Locale;

import org.junit.Assert;
import org.junit.Test;

import drafts.AndersenBroadieUpperBoundEstimation;
import drafts.BermudanSwaption;
import drafts.BermudanSwaptionValueEstimatorInterface;
import drafts.DeltaHedgingUpperBound;
import drafts.SimpleLowerBoundEstimation;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;

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
	static double tolerance = 0.01; // should be tightened pending further improvement
	static int numberOfExercisePeriods = 10;
	private static int numberOfSubsimulationsStepA=1000;
	private static int numberOfSubsimulationsStepB=1000;
	

	@Test
	public static void testSwaptionValuationMethod() throws CalculationException {

		LIBORModelMonteCarloSimulationModel liborModel = CreateTestModel.createLIBORMarketModel();
		// print head of comparison table
		System.out.println("Bermudan Swaption prices:\n");
		System.out.println(
				"FirstFixingDate\tLower Bound\tUpper Bound(AB)\tUpperBound(Deltas)\tDeviation(AB)\tDeviation(Delta subsimfree)");
		// "EvaluationDate Lower Bound Upper Bound(AB) Deviation(AB) ");


		for (int maturityIndex = 1; maturityIndex < liborModel.getNumberOfLibors()
				- numberOfExercisePeriods; maturityIndex++) {
			double firstFixingDate = liborModel.getLiborPeriod(maturityIndex);
			CreateTestBermudanSwaption.setFirstFixingDate(firstFixingDate);
			CreateTestBermudanSwaption.setNumberOfExercisePeriods(numberOfExercisePeriods);
			BermudanSwaption testSwaption = CreateTestBermudanSwaption.createBermudanSwaption();

			/*
			 * Value a bermudan swaption
			 */

			System.out.print(formatterTime.format(firstFixingDate) + "\t");
			

			SimpleLowerBoundEstimation lowerBound = new SimpleLowerBoundEstimation();

			double lowerBoundValue = timingValuationTest(lowerBound, testSwaption, liborModel);

			// AB upper bound approximation
			AndersenBroadieUpperBoundEstimation ABupperBound = new AndersenBroadieUpperBoundEstimation(lowerBound,
					numberOfSubsimulationsStepA, numberOfSubsimulationsStepB);
			double ABupperBoundValue = timingValuationTest(ABupperBound, testSwaption, liborModel);
			// Upper bound using martingale construction via Delta hedging
			DeltaHedgingUpperBound DeltaUpperBound = new DeltaHedgingUpperBound(lowerBound);
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
	
	//help method for calculation and printing of values and calculation time
	private static double timingValuationTest (BermudanSwaptionValueEstimatorInterface valuationMethod, BermudanSwaption swaption, LIBORModelMonteCarloSimulationModel liborModel) throws CalculationException
	{
		swaption.setValuationMethod(valuationMethod);
		long timingValuationStart = System.currentTimeMillis();
		double simulatedValue = swaption.getValue(liborModel);
		long timingValuationEnd = System.currentTimeMillis();
		double lastOperationTimingValuation = (timingValuationEnd - timingValuationStart) / 1000.0;
		System.out.print(formatterValue.format(simulatedValue)+"\t");
		System.out.print(formattterTime.format(lastOperationTimingValuation)+ "\t");
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