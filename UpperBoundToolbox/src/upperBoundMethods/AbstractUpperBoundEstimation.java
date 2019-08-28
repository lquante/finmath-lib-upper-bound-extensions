package upperBoundMethods;

import java.util.ArrayList;

import bermudanSwaptionFramework.BermudanSwaption;
import bermudanSwaptionFramework.BermudanSwaptionValueEstimatorInterface;
import lowerBoundMethods.AbstractLowerBoundEstimationInputForUpperBound;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.conditionalexpectation.RegressionBasisFunctionsProvider;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;

/**
 * This class provides the framework for an upper bound estimation based on an
 * abstract estimation of a martingale process following the primal-dual
 * approach.
 * 
 * The abstract method calculateMartingaleApproximation needs to be implemented
 * depending on the specific method.
 * 
 * @author Lennart Quante
 * @version 1.0
 *
 * 
 */
public abstract class AbstractUpperBoundEstimation implements BermudanSwaptionValueEstimatorInterface {
	// basis functions used in the lower bound
	protected RegressionBasisFunctionsProvider regressionBasisFunctionsProvider;
	// product to be evaluated
	BermudanSwaption bermudanSwaption;
	private AbstractLowerBoundEstimationInputForUpperBound lowerBoundMethod;

	private double weightOfMartingale;

	/**
	 * Constructor setting basic properties of the upper bound estimation.
	 *@param lowerBoundMethod The lower bound method to be used as input for the
	 *                         upper bound approximation.
	 *@param weightOfMartingale The weight with which the martingale approximation should be added. 1=upper bound, 0 results in lower bound, 0.5= point wise estimate of A-B.                  
	 */
	public AbstractUpperBoundEstimation(AbstractLowerBoundEstimationInputForUpperBound lowerBoundMethod, double weightOfMartingale) {

		this.lowerBoundMethod = lowerBoundMethod;
		this.weightOfMartingale = weightOfMartingale;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * drafts.BermudanSwaptionValueEstimatorInterface#getValueEstimation(drafts.
	 * BermudanSwaption, double,
	 * net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel)
	 */
	@Override
	public RandomVariable getValueEstimation(BermudanSwaption bermudanOption, double evaluationTime,
			LIBORModelMonteCarloSimulationModel model, RandomVariable triggerValuesInput) throws CalculationException {

		this.bermudanSwaption = bermudanOption;
		RandomVariable optionValue = model.getRandomVariableForConstant(0);
		// calculate value of lower bound

		this.bermudanSwaption.setValuationMethod(lowerBoundMethod);
		this.bermudanSwaption.getValue(evaluationTime, model);
		RandomVariable[] cacheUnderlying = lowerBoundMethod.getCacheValuesOfUnderlying();
		RandomVariable[] cacheOptionValues = lowerBoundMethod.getCacheOptionValues();
		// Check exercise condition of lower bound method
		RandomVariable[] cacheTriggers = lowerBoundMethod.getCacheTriggerValues();
		// calculate option time index for evaluation
		int numberOfOptionPeriods = this.bermudanSwaption.getFixingDates().length-1;
		double firstFixingDate = this.bermudanSwaption.getFixingDates()[0];
		int optionTimeIndex;
		if (evaluationTime <= firstFixingDate)
			optionTimeIndex = 0;
		else
			for (optionTimeIndex = 1; evaluationTime > this.bermudanSwaption
					.getFixingDates()[optionTimeIndex]; optionTimeIndex++)
				;
		// calculate martingale approximation

		ArrayList <RandomVariable> martingaleArrays = calculateMartingaleApproximation(optionTimeIndex, model, cacheUnderlying,
				cacheOptionValues, cacheTriggers);

		RandomVariable martingaleApproximation = model.getRandomVariableForConstant(0);

		
		// calculate the maximum of the estimated martingale for each remaining period
		for (int optionPeriod = optionTimeIndex; optionPeriod < numberOfOptionPeriods; optionPeriod++)
			martingaleApproximation = martingaleApproximation.floor(cacheUnderlying[optionPeriod].sub(martingaleArrays.get(optionPeriod)));
		// Note that values is a relative price - no numeraire division is required
		RandomVariable numeraireAtEvaluationTime = model.getNumeraire(evaluationTime);
		RandomVariable monteCarloProbabilitiesAtEvaluationTime = model.getMonteCarloWeights(evaluationTime);
		RandomVariable discountFactor = numeraireAtEvaluationTime.div(monteCarloProbabilitiesAtEvaluationTime);
		optionValue = cacheOptionValues[optionTimeIndex].mult(discountFactor).add(martingaleApproximation).mult(weightOfMartingale);
		return optionValue;

	}
	
	/**
	 * This method needs to be specified to estimate the martingale component of the
	 * upper bound estimation using the following parameters:
	 * 
	 * @param period            period of the evaluation time in the model time
	 *                          discretization
	 * @param model             the LIBORModelMonteCarloSimulationModel to be used
	 * @param cacheUnderlying   The cached exercise values from the lower bound
	 *                          method.
	 * @param cacheOptionValues The cached option values from the lower bound
	 *                          method.
	 * @param triggerValues     The cached trigger values from the lower bound
	 *                          method.
	 * @return The Monte Carlo estimate of the martingale component, as a double.
	 * @throws CalculationException
	 */
	protected abstract ArrayList<RandomVariable> calculateMartingaleApproximation(int period, LIBORModelMonteCarloSimulationModel model,
			RandomVariable[] cacheUnderlying, RandomVariable[] cacheOptionValues, RandomVariable[] triggerValues)
					throws CalculationException;




	// some getters and setters

	/**
	 * @return the weightOfMartingale
	 */
	public double getWeightOfMartingale() {
		return weightOfMartingale;
	}

	/**
	 * @param weightOfMartingale the weightOfMartingale to set
	 */
	public void setWeightOfMartingale(double weightOfMartingale) {
		this.weightOfMartingale = weightOfMartingale;
	}

}