/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 28 Feb 2015
 */
package tests;

import org.junit.Test;

import bermudanswaptionframework.BermudanSwaption;
import bermudanswaptionframework.BermudanSwaptionValueEstimatorInterface;
import org.junit.Assert;
import lowerboundmethods.SimpleLowerBoundEstimation;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import upperboundmethods.AndersenBroadieUpperBoundEstimation;
import upperboundmethods.DeltaHedgingUpperBound;

/**
 * This class tests if a valuation of a swaption using the different valuation methods suceeds:
 *
 * @author Lennart Quante
 */

public class TestBermudanSwaption {
	
	
	
	private int numberOfExercisePeriods;
	private double periodLength;
	private double swaprate;

	

	@Test
	public void testLowerBoundSwaptionValuation() throws CalculationException {
		
		// set parameters
		BermudanSwaptionValueEstimatorInterface valuationMethod = new SimpleLowerBoundEstimation();
		double evaluationTime = 2;
		double swaptionValue = evaluateSwaption(evaluationTime,valuationMethod);
		System.out.println("Lower Bound Bermudan Swaption value at time " + evaluationTime + ": " + "" + swaptionValue);
		Assert.assertNotNull(swaptionValue);
	}
	
	@Test
	public void testAndersenBroadieSwaptionValuation() throws CalculationException {
		
		// set parameters
		BermudanSwaptionValueEstimatorInterface valuationMethod = new AndersenBroadieUpperBoundEstimation(100);
		double evaluationTime = 2;
		double swaptionValue = evaluateSwaption(evaluationTime,valuationMethod);
		System.out.println("AB Bermudan Swaption value at time " + evaluationTime + ": " + "" + swaptionValue);
		Assert.assertNotNull(swaptionValue);
	}
	
	@Test
	public void testDeltaUpperBoundSwaptionValuation() throws CalculationException {
		
		// set parameters
		BermudanSwaptionValueEstimatorInterface valuationMethod = new DeltaHedgingUpperBound();
		double evaluationTime = 2;
		double swaptionValue = evaluateSwaption(evaluationTime,valuationMethod);
		System.out.println("Delta Upper Bound Bermudan Swaption value at time " + evaluationTime + ": " + "" + swaptionValue);
		Assert.assertNotNull(swaptionValue);
	}
	
	
	private double evaluateSwaption (double evaluationTime, BermudanSwaptionValueEstimatorInterface valuationMethod) throws CalculationException
	{
		// create swaption
				double firstFixingDate = 1;
				this.numberOfExercisePeriods=10;
				this.periodLength =1;
				this.swaprate = 0.02;
				BermudanSwaption testSwaption = constructBermudanSwaption(firstFixingDate, valuationMethod);
				// create model as specified in model test class:
				LIBORModelMonteCarloSimulationModel testModel = CreateTestModel.createLIBORMarketModel();
				// get value
				return testSwaption.getValue(evaluationTime, testModel).getAverage();
	}
	
	
	
	/**
	 * Method to construct callable Bermudan swaption for testing purposes from reduced input. All fixing dates are exercise dates, equal period length, notional ==1 and constant swaprate
	 * @param firstFixingDate first fixing date of the swaption
	 * @param numberOfExercisePeriods
	 * @param periodLength
	 * @param swaprate the constant swaprate
	 * @param valuationMethod valuation method to be used
	 * @return A Bermudan Swaption with the given parameters
	 */
	public BermudanSwaption constructBermudanSwaption (double firstFixingDate,BermudanSwaptionValueEstimatorInterface valuationMethod){
		
		double[] fixingDates = new double[numberOfExercisePeriods];
		double[] periodLengths = new double[numberOfExercisePeriods];
		double[] paymentDates = new double[numberOfExercisePeriods];
		double[] notionals = new double[numberOfExercisePeriods + 1];
		boolean[] isPeriodStartDateExerciseDate = new boolean[numberOfExercisePeriods + 1];
		double[] swapTenor = new double[numberOfExercisePeriods]; // to be passed to the analytical approximation method

		for (int periodStartIndex = 0; periodStartIndex < numberOfExercisePeriods; periodStartIndex++) {
			fixingDates[periodStartIndex] = firstFixingDate + periodStartIndex * periodLength;
			periodLengths[periodStartIndex] = periodLength;
			paymentDates[periodStartIndex] = firstFixingDate + (periodStartIndex + 1) * periodLength;
			swapTenor[periodStartIndex] = periodLength;
			notionals[periodStartIndex] = 1;
			isPeriodStartDateExerciseDate[periodStartIndex] = true;
			
		}
		notionals[numberOfExercisePeriods] = 1;
		isPeriodStartDateExerciseDate[numberOfExercisePeriods] = true;

		// Set swap rates for each period
		double[] swaprates = new double[numberOfExercisePeriods];
		for (int periodStartIndex = 0; periodStartIndex < numberOfExercisePeriods; periodStartIndex++) {
			swaprates[periodStartIndex] = swaprate;
		}
		return new BermudanSwaption("EURO", isPeriodStartDateExerciseDate, fixingDates, periodLengths, paymentDates, notionals, true,swaprates,valuationMethod);
	
	}

	/**
	 * @return the numberOfExercisePeriods
	 */
	public int getNumberOfExercisePeriods() {
		return numberOfExercisePeriods;
	}

	/**
	 * @param numberOfExercisePeriods the numberOfExercisePeriods to set
	 */
	public void setNumberOfExercisePeriods(int numberOfExercisePeriods) {
		this.numberOfExercisePeriods = numberOfExercisePeriods;
	}

	/**
	 * @return the periodLength
	 */
	public double getPeriodLength() {
		return periodLength;
	}

	/**
	 * @param periodLength the periodLength to set
	 */
	public void setPeriodLength(double periodLength) {
		this.periodLength = periodLength;
	}

	/**
	 * @return the swaprate
	 */
	public double getSwaprate() {
		return swaprate;
	}

	/**
	 * @param swaprate the swaprate to set
	 */
	public void setSwaprate(double swaprate) {
		this.swaprate = swaprate;
	}

}
