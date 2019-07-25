/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 28 Feb 2015
 */
package tests;

import org.junit.Test;

import drafts.BermudanSwaption;
import drafts.BermudanSwaptionValueEstimatorInterface;
import drafts.SimpleLowerBoundEstimation;
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
	static int numberOfPeriods = 10; // number of (possible) exercise dates
	static double swapPeriodLength = 0.5;
	static String currency = "EURO";
	static double swaprate = 0.02; // getParSwaprate(liborModel, swapTenor);
	static BermudanSwaption bermudanSwaption;
	static BermudanSwaptionValueEstimatorInterface valuationInterface = new SimpleLowerBoundEstimation(); // valuation
	// interface
	static double firstFixingDate = 0.5;

	static double evaluationTime = 0;
	// model field
	LIBORModelMonteCarloSimulationModel model;

	public CreateTestBermudanSwaption() throws CalculationException {
		model = CreateTestModel.createLIBORMarketModel();
		this.bermudanSwaption = createBermudanSwaption();
	}

	@Test
	public void testSwaptionValuation() throws CalculationException {
		// test valuation
		double swaptionValue = bermudanSwaption.getValue(evaluationTime, model).getAverage();
		System.out.println("Bermudan Swaption value at time " + evaluationTime + ": " + "" + swaptionValue);

	}

	public static BermudanSwaption createBermudanSwaption() {
		// Create a rudimental bermudan swaption

		double[] fixingDates = new double[numberOfPeriods];
		double[] paymentDates = new double[numberOfPeriods];

		double[] notionals = new double[numberOfPeriods + 1];
		boolean[] startingDate = new boolean[numberOfPeriods + 1];
		double[] swapTenor = new double[numberOfPeriods]; // to be passed to the analytical approximation method

		for (int periodStartIndex = 0; periodStartIndex < numberOfPeriods; periodStartIndex++) {
			fixingDates[periodStartIndex] = firstFixingDate + periodStartIndex * swapPeriodLength;
			paymentDates[periodStartIndex] = firstFixingDate + (periodStartIndex + 1) * swapPeriodLength;
			swapTenor[periodStartIndex] = swapPeriodLength;
			notionals[periodStartIndex] = 1;
			startingDate[periodStartIndex] = true;
		}
		notionals[numberOfPeriods] = 1;
		startingDate[numberOfPeriods] = true;

		// Set swap rates for each period
		double[] swaprates = new double[numberOfPeriods];
		for (int periodStartIndex = 0; periodStartIndex < numberOfPeriods; periodStartIndex++) {
			swaprates[periodStartIndex] = swaprate;
		}

		return new BermudanSwaption(currency, startingDate, fixingDates, swapTenor, paymentDates, notionals, true,
				swaprates, valuationInterface);
	}

	// some getters and setters

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

}
