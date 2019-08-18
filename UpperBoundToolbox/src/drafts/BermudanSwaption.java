package drafts;

import java.util.Arrays;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;

public class BermudanSwaption extends AbstractLIBORBermudanOption {

	private final double[] swaprates; // Vector of strikes
	private BermudanSwaptionValueEstimatorInterface valuationMethod;

	/**
	 * @param swaprates       Array of swaprates
	 * @param valuationMethod The valuation method to be used.
	 */
	public BermudanSwaption(String currency, boolean[] isPeriodStartDateExerciseDate, double[] fixingDates,
			double[] periodLengths, double[] paymentDates, double[] periodNotionals, boolean isCallable,
			double[] swaprates, BermudanSwaptionValueEstimatorInterface valuationMethod) {
		super(currency, isPeriodStartDateExerciseDate, fixingDates, periodLengths, paymentDates, periodNotionals,
				isCallable);
		this.swaprates = swaprates;
		this.valuationMethod = valuationMethod;
	}

	public BermudanSwaption getBermudanSwaptionWithChangedValuationMethod(
			BermudanSwaptionValueEstimatorInterface valuationMethod) {
		return new BermudanSwaption(this.getCurrency(), this.getIsPeriodStartDateExerciseDate(), this.getFixingDates(),
				this.getPeriodLengths(), this.getPaymentDates(), this.getPeriodNotionals(), this.isCallable(),
				this.swaprates, valuationMethod);
	}

	@Override
	public RandomVariable getValue(double evaluationTime, LIBORModelMonteCarloSimulationModel model)
			throws CalculationException {
		return valuationMethod.getValueEstimation(this, evaluationTime, model, null);
	}

	// method to produce clone with same properties, but later starting date:

	/**
	 * @param startingIndex  only from the exisiting fixing dates)
	 * @return
	 */
	public BermudanSwaption getCloneWithModifiedStartingPeriod(int startingIndex) {
		
		if (startingIndex >=this.getFixingDates().length)
			startingIndex=this.getFixingDates().length-1;
		
		// shorten all arrays for a new BermudanSwaption
		
		double[] adjustedFixingDates = Arrays.copyOfRange(this.getFixingDates(), startingIndex,
				this.getFixingDates().length);
		boolean[] adjustedIsPeriodStartDateExerciseDate = Arrays.copyOfRange(this.getIsPeriodStartDateExerciseDate(),
				startingIndex, this.getFixingDates().length);
		double[] adjustedPeriodLengths = Arrays.copyOfRange(this.getPeriodLengths(), startingIndex,
				this.getFixingDates().length);
		double[] adjustedPaymentDates = Arrays.copyOfRange(this.getPaymentDates(), startingIndex,
				this.getFixingDates().length);
		double[] adjustedPeriodNotionals = Arrays.copyOfRange(this.getPeriodNotionals(), startingIndex,
				this.getFixingDates().length);
		double[] adjustedSwapRates = Arrays.copyOfRange(this.getSwaprates(), startingIndex,
				this.getFixingDates().length);
		BermudanSwaptionValueEstimatorInterface unadjustedValuationMethod = this.getValuationMethod();
		return new BermudanSwaption(this.getCurrency(), adjustedIsPeriodStartDateExerciseDate, adjustedFixingDates,
				adjustedPeriodLengths, adjustedPaymentDates, adjustedPeriodNotionals, this.isCallable(),
				adjustedSwapRates, unadjustedValuationMethod);

	}

	public BermudanSwaption getCloneWithModifiedStartingAndFinalPeriod(int startingIndex,int finalIndex) {
		
		if (startingIndex >=this.getFixingDates().length)
			startingIndex=this.getFixingDates().length-1;
		finalIndex+=1;
		if (finalIndex >=this.getFixingDates().length)
			finalIndex=this.getFixingDates().length;
		
		// shorten all arrays for a new BermudanSwaption

		double[] adjustedFixingDates = Arrays.copyOfRange(this.getFixingDates(), startingIndex, finalIndex);
		boolean[] adjustedIsPeriodStartDateExerciseDate = Arrays.copyOfRange(this.getIsPeriodStartDateExerciseDate(), startingIndex,
				finalIndex);
		double[] adjustedPeriodLengths = Arrays.copyOfRange(this.getPeriodLengths(), startingIndex, finalIndex);
		double[] adjustedPaymentDates = Arrays.copyOfRange(this.getPaymentDates(), startingIndex, finalIndex);
		double[] adjustedPeriodNotionals = Arrays.copyOfRange(this.getPeriodNotionals(), startingIndex, finalIndex);
		double[] adjustedSwapRates = Arrays.copyOfRange(this.getSwaprates(), startingIndex, finalIndex);
		
		return new BermudanSwaption(this.getCurrency(), adjustedIsPeriodStartDateExerciseDate, adjustedFixingDates,
				adjustedPeriodLengths, adjustedPaymentDates, adjustedPeriodNotionals, this.isCallable(),
				adjustedSwapRates, this.getValuationMethod());

	}

	// some getters

	public BermudanSwaptionValueEstimatorInterface getValuationMethod() {
		return valuationMethod;
	}

	public void setValuationMethod(BermudanSwaptionValueEstimatorInterface valuationMethod) {
		this.valuationMethod = valuationMethod;
	}

	public double[] getSwaprates() {
		return swaprates;

	}
}
