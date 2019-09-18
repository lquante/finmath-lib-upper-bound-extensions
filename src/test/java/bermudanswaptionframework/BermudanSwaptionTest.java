/**
 * 
 */
package bermudanswaptionframework;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import bermudanswaptionframework.BermudanSwaption;
import org.junit.Assert;
import lowerboundmethods.SimpleLowerBoundEstimation;
import upperboundmethods.DeltaHedgingUpperBound;

/**
 * @author Lennart
 *
 */
public class BermudanSwaptionTest {

	
	// parameters to be used:
	int numberOfPeriods = 5;
	
	
	/**
	 * Test method for {@link bermudanswaptionframework.BermudanSwaption#getBermudanSwaptionWithChangedValuationMethod(bermudanswaptionframework.BermudanSwaptionValueEstimatorInterface)}.
	 */
	@Test
	public void testGetBermudanSwaptionWithChangedValuationMethod() {
		BermudanSwaption testSwaption = constructTestSwaption();
		BermudanSwaptionValueEstimatorInterface anotherValuationMethod = new DeltaHedgingUpperBound();
		BermudanSwaption changedSwaption = testSwaption.getBermudanSwaptionWithChangedValuationMethod(anotherValuationMethod);
		Assert.assertTrue(changedSwaption.getValuationMethod()==anotherValuationMethod);
	}

	/**
	 * Test method for {@link bermudanswaptionframework.BermudanSwaption#getCloneWithModifiedStartingPeriod(int)}.
	 */
	@Test
	public void testGetCloneWithModifiedStartingPeriod() {
		BermudanSwaption testSwaption = constructTestSwaption();
		int startingPeriod = numberOfPeriods-numberOfPeriods/2;
		BermudanSwaption changedSwaption = testSwaption.getCloneWithModifiedStartingPeriod(startingPeriod);
		Assert.assertTrue(changedSwaption.getFixingDates()[0]==startingPeriod);
	}

	/**
	 * Test method for {@link bermudanswaptionframework.BermudanSwaption#getCloneWithModifiedStartingAndFinalPeriod(int, int)}.
	 */
	@Test
	public void testGetCloneWithModifiedStartingAndFinalPeriod() {
		BermudanSwaption testSwaption = constructTestSwaption();
		int startingPeriod = numberOfPeriods-numberOfPeriods/2;
		int finalPeriod = numberOfPeriods-1;
		BermudanSwaption changedSwaption = testSwaption.getCloneWithModifiedStartingAndFinalPeriod(startingPeriod, finalPeriod);
		double[] fixingDates = changedSwaption.getFixingDates();
		int numberOfNewPeriods = fixingDates.length;
		Assert.assertTrue(fixingDates[0]==startingPeriod & fixingDates[numberOfNewPeriods-1]==finalPeriod );
	}
	
	
	/** help method to create swaption
	 * @return a Bermudan swaption for testing with numberOfPeriods periods
	*/
	BermudanSwaption constructTestSwaption () {
		boolean isCallable = true;
		boolean[] isPeriodStartDateExerciseDate = new boolean[numberOfPeriods];
		double[] fixingDates = new double [numberOfPeriods];
		double[] periodLengths = new double [numberOfPeriods];
		double[] paymentDates = new double [numberOfPeriods];
		double[] periodNotionals =  new double [numberOfPeriods];
		double[] swaprates = new double [numberOfPeriods];
		for(int i=0; i<numberOfPeriods;i++)
		{
			isPeriodStartDateExerciseDate[i]=true;
			fixingDates[i]= i;	
			periodLengths[i]= 1;	
			paymentDates[i] = fixingDates[i]+periodLengths[i];
			periodNotionals[i] = 1;
			swaprates[i]=0.02;
		}
		BermudanSwaptionValueEstimatorInterface valuationMethod = new SimpleLowerBoundEstimation();
		return new BermudanSwaption(isPeriodStartDateExerciseDate,fixingDates,periodLengths,paymentDates,periodNotionals, isCallable, swaprates, valuationMethod);
	}
	

}
