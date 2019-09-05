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
import simulationMethods.TestModelFactory;
import simulationMethods.TestSwaptionFactory;
import upperboundmethods.AndersenBroadieUpperBoundEstimation;
import upperboundmethods.DeltaHedgingUpperBound;

/**
 * This class tests if a valuation of a swaption using the different valuation methods suceeds:
 *
 * @author Lennart Quante
 */

public class TestBermudanSwaptionValuation {
	
	
	
	private int numberOfExercisePeriods=10;
	private double periodLength=1;
	private double swaprate=0.02;

	

	@Test
	public void testLowerBoundSwaptionValuation() throws CalculationException {
		
		// set parameters
		BermudanSwaptionValueEstimatorInterface valuationMethod = new SimpleLowerBoundEstimation();
		double evaluationTime = 1;
		double swaptionValue = evaluateSwaption(evaluationTime,valuationMethod);
		System.out.println("Lower Bound Bermudan Swaption value at time " + evaluationTime + ": " + "" + swaptionValue);
		Assert.assertNotNull(swaptionValue);
	}
	
	@Test
	public void testAndersenBroadieSwaptionValuation() throws CalculationException {
		
		// set parameters
		BermudanSwaptionValueEstimatorInterface valuationMethod = new AndersenBroadieUpperBoundEstimation(100);
		double evaluationTime = 1;
		double swaptionValue = evaluateSwaption(evaluationTime,valuationMethod);
		System.out.println("AB Bermudan Swaption value at time " + evaluationTime + ": " + "" + swaptionValue);
		Assert.assertNotNull(swaptionValue);
	}
	
	@Test
	public void testDeltaUpperBoundSwaptionValuation() throws CalculationException {
		
		// set parameters
		BermudanSwaptionValueEstimatorInterface valuationMethod = new DeltaHedgingUpperBound();
		double evaluationTime = 1;
		double swaptionValue = evaluateSwaption(evaluationTime,valuationMethod);
		System.out.println("Delta Upper Bound Bermudan Swaption value at time " + evaluationTime + ": " + "" + swaptionValue);
		Assert.assertNotNull(swaptionValue);
	}
	
	
	private double evaluateSwaption (double evaluationTime, BermudanSwaptionValueEstimatorInterface valuationMethod) throws CalculationException
	{
		// create swaption "factory"
		TestSwaptionFactory testSwaptionProvider = new TestSwaptionFactory(numberOfExercisePeriods,periodLength,swaprate);
		
				double firstFixingDate = 1;
				this.numberOfExercisePeriods=10;
				this.periodLength =1;
				this.swaprate = 0.02;
				BermudanSwaption testSwaption = testSwaptionProvider.constructBermudanSwaption(firstFixingDate, valuationMethod);
				// create model as specified in model test class:
				LIBORModelMonteCarloSimulationModel testModel = TestModelFactory.createLIBORMarketModel();
				// get value
				return testSwaption.getValue(evaluationTime, testModel).getAverage();
	}
}
	
	
	
