package tests;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

import org.junit.Assert;
import org.junit.Test;

import drafts.AndersenBroadieUpperBoundEstimation;
import drafts.BermudanSwaption;
import drafts.DeltaHedgingUpperBound;
import drafts.SimpleLowerBoundEstimation;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;

public class TestValuationMethods {
	// formatter for values

	private static DecimalFormat formatterTime = new DecimalFormat("00.00", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterValue = new DecimalFormat(" ##0.00000%;-##0.00000%",
			new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterDeviation = new DecimalFormat(" 0.00000E00;-0.00000E00",
			new DecimalFormatSymbols(Locale.ENGLISH));

	// set tolerance for difference between upper and lower bound methods
	double tolerance = 0.01; // should be tightened pending further improvement
	int numberOfPeriods = 10;
	private int numberOfSubsimulationsStepA=100;
	private int numberOfSubsimulationsStepB=100;

	@Test
	public void testSwaptionValuationMethod() throws CalculationException {

		LIBORModelMonteCarloSimulationModel liborModel = CreateTestModel.createLIBORMarketModel();
		// print head of comparison table
		System.out.println("Bermudan Swaption prices:\n");
		System.out.println(
				"FirstFixingDate Lower Bound \t Upper Bound(AB) \t UpperBound(Deltas)\tDeviation(AB)\tDeviation(Delta subsimfree)");
		// "EvaluationDate Lower Bound Upper Bound(AB) Deviation(AB) ");


		for (int maturityIndex = 1; maturityIndex < liborModel.getNumberOfLibors()
				- numberOfPeriods; maturityIndex++) {
			double firstFixingDate = liborModel.getLiborPeriod(maturityIndex);

			CreateTestBermudanSwaption.setFirstFixingDate(firstFixingDate);
			BermudanSwaption testSwaption = CreateTestBermudanSwaption.createBermudanSwaption();

			/*
			 * Value a bermudan swaption
			 */

			System.out.print(formatterTime.format(firstFixingDate) + "\t\t");
			

			SimpleLowerBoundEstimation lowerBound = new SimpleLowerBoundEstimation();

			testSwaption.setValuationMethod(lowerBound);
			double lowerBoundValue = testSwaption.getValue(liborModel);
			System.out.print(formatterValue.format(lowerBoundValue) + "          ");

			// AB upper bound approximation
			AndersenBroadieUpperBoundEstimation ABupperBound = new AndersenBroadieUpperBoundEstimation(lowerBound,
					numberOfSubsimulationsStepA, numberOfSubsimulationsStepB);
			testSwaption.setValuationMethod(ABupperBound);
			double ABupperBoundValue = testSwaption.getValue(liborModel);
			System.out.print(formatterValue.format(ABupperBoundValue) + "\t");

			// Upper bound using martingale construction via Delta hedging
			DeltaHedgingUpperBound DeltaUpperBound = new DeltaHedgingUpperBound(lowerBound);
			testSwaption.setValuationMethod(DeltaUpperBound);
			double DeltaUpperBoundValue = testSwaption.getValue(liborModel);
			System.out.print(formatterValue.format(DeltaUpperBoundValue) + "\t");

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

}
