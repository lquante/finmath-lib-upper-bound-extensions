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
	int numberOfPeriods = 10; // number of (possible) exercise dates
	double swapPeriodLength = 0.5;
	String currency = "EURO";
	double swaprate = 0.02; // getParSwaprate(liborModel, swapTenor);
	private BermudanSwaption bermudanSwaption;
	BermudanSwaptionValueEstimatorInterface valuationInterface = new SimpleLowerBoundEstimation(); // valuation
																									// interface
	double firstFixingDate = 0.5;
	double evaluationTime = 0;
	// model field
	LIBORModelMonteCarloSimulationModel model;

	public CreateTestBermudanSwaption() throws CalculationException {
		model = CreateTestModel.createLIBORMarketModel();
		// Create a rudimental bermudan swaption

		double[] fixingDates = new double[numberOfPeriods];
		double[] paymentDates = new double[numberOfPeriods]; // simply a merge of fixing and payment dates (which
																// obviously overlap)

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

		this.bermudanSwaption = new BermudanSwaption(currency, startingDate, fixingDates, swapTenor, paymentDates,
				notionals, true, swaprates, valuationInterface);
	}

	@Test
	public void testSwaptionValuation() throws CalculationException {
		// test valuation
		double swaptionValue = bermudanSwaption.getValue(evaluationTime, model).getAverage();
		System.out.println("Bermudan Swaption value at time " + evaluationTime + ": " + "" + swaptionValue);

	}
}
