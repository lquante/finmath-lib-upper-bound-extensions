/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 28 Feb 2015
 */
package tests;

import org.junit.Test;

import bermudanswaptionframework.BermudanSwaption;
import bermudanswaptionframework.BermudanSwaptionValueEstimatorInterface;
import lowerboundmethods.SimpleLowerBoundEstimation;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;

/**
 * This class tests the LIBOR market model and products.
 *
 * @author Christian Fries
 * @author Lorenzo Torricelli
 */

public class CreateTestBermudanSwaption {
	// option parameters
	static int numberOfExercisePeriods = 10; // number of (possible) exercise dates
	static double swapPeriodLength = 1;
	static String currency = "EURO";
	static double swaprate = 0.021; // getParSwaprate(liborModel, swapTenor);
	static BermudanSwaption bermudanSwaption;
	static BermudanSwaptionValueEstimatorInterface valuationInterface = new SimpleLowerBoundEstimation(); // valuation
	// interface
	static double firstFixingDate = 0.5;

	static double evaluationTime = 0;
	// model field
	LIBORModelMonteCarloSimulationModel model;

	public CreateTestBermudanSwaption() throws CalculationException {
		model = CreateTestModel.createLIBORMarketModel();
		CreateTestBermudanSwaption.bermudanSwaption = createBermudanSwaption();
	}

	@Test
	public void testSwaptionValuation() throws CalculationException {
		// test valuation
		double swaptionValue = bermudanSwaption.getValue(evaluationTime, model).getAverage();
		System.out.println("Bermudan Swaption value at time " + evaluationTime + ": " + "" + swaptionValue);

	}

	public static BermudanSwaption createBermudanSwaption() {
		// Create a rudimental bermudan swaption

		double[] fixingDates = new double[numberOfExercisePeriods];
		double[] paymentDates = new double[numberOfExercisePeriods];

		double[] notionals = new double[numberOfExercisePeriods + 1];
		boolean[] startingDate = new boolean[numberOfExercisePeriods + 1];
		double[] swapTenor = new double[numberOfExercisePeriods]; // to be passed to the analytical approximation method

		for (int periodStartIndex = 0; periodStartIndex < numberOfExercisePeriods; periodStartIndex++) {
			fixingDates[periodStartIndex] = firstFixingDate + periodStartIndex * swapPeriodLength;
			paymentDates[periodStartIndex] = firstFixingDate + (periodStartIndex + 1) * swapPeriodLength;
			swapTenor[periodStartIndex] = swapPeriodLength;
			notionals[periodStartIndex] = 1;
			startingDate[periodStartIndex] = true;
		}
		notionals[numberOfExercisePeriods] = 1;
		startingDate[numberOfExercisePeriods] = true;

		// Set swap rates for each period
		double[] swaprates = new double[numberOfExercisePeriods];
		for (int periodStartIndex = 0; periodStartIndex < numberOfExercisePeriods; periodStartIndex++) {
			swaprates[periodStartIndex] = swaprate;
		}

		return new BermudanSwaption(currency, startingDate, fixingDates, swapTenor, paymentDates, notionals, true,
				swaprates, valuationInterface);
	}

	// some getters and setters

	/**
	 * @param numberOfPeriods the numberOfPeriods to set
	 */
	public static void setNumberOfExercisePeriods(int numberOfPeriods) {
		CreateTestBermudanSwaption.numberOfExercisePeriods = numberOfPeriods;
	}

	/**
	 * @return firstFixingDate
	 */
	public static double getFirstFixingDate() {
		return firstFixingDate;
	}

	/**
	 * @param firstFixingDate the firstFixingDate to set
	 */
	public static void setFirstFixingDate(double firstFixingDate) {
		CreateTestBermudanSwaption.firstFixingDate = firstFixingDate;
	}

	public static void setPeriodLength(double periodLength) {
		
		CreateTestBermudanSwaption.swapPeriodLength = periodLength;
	}

	/**
	 * @return the swaprate
	 */
	public static double getSwaprate() {
		return swaprate;
	}

	/**
	 * @param swaprate the swaprate to set
	 */
	public static void setSwaprate(double swaprate) {
		CreateTestBermudanSwaption.swaprate = swaprate;
	}

}
