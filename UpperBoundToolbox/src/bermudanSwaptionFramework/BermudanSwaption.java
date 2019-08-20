package bermudanSwaptionFramework;

import java.util.Arrays;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;

/**
 * Implements a Bermudan swaption product..
 * @author Lennart Quante
 * @version 1.0
 */
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

	/**
	 * @param valuationMethod
	 * @return A Bermudan Swaption with the input valuation method
	 */
	public BermudanSwaption getBermudanSwaptionWithChangedValuationMethod(
			BermudanSwaptionValueEstimatorInterface valuationMethod) {
		return new BermudanSwaption(this.getCurrency(), this.getIsPeriodStartDateExerciseDate(), this.getFixingDates(),
				this.getPeriodLengths(), this.getPaymentDates(), this.getPeriodNotionals(), this.isCallable(),
				this.swaprates, valuationMethod);
	}

	/**
	 * Method using the valuation method specified in the option.
	 */
	@Override
	public RandomVariable getValue(double evaluationTime, LIBORModelMonteCarloSimulationModel model)
			throws CalculationException {
		return valuationMethod.getValueEstimation(this, evaluationTime, model, null);
	}

	// method to produce clone with same properties, but later starting date:

	/**
	 * Shifts the first exercise date of an Bermudan swaption needed e.g. in method of Andersen-Broadie.
	 * @param startingIndex  first fixing date to be kept
	 * @return Clone of the BermudanSwaption with the input starting index.
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

	/**
	 *   Shifts the first and final exercise date of an Bermudan swaption needed e.g. in method of Andersen-Broadie.
	 * @param startingIndex  first fixing date to be kept
	 * @param finalIndex last fixing date to be kept (inclusive)
	 * @return Clone of the BermudanSwaption with the input starting index and final incex.
	 */
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

	/**
	 * @return the valuation method of the Bermudan swaption
	 */
	public BermudanSwaptionValueEstimatorInterface getValuationMethod() {
		return valuationMethod;
	}

	/**
	 * Changes the valuation method of an Bermudan swaption without generating a new object.
	 * @param valuationMethod
	 */
	public void setValuationMethod(BermudanSwaptionValueEstimatorInterface valuationMethod) {
		this.valuationMethod = valuationMethod;
	}

	/**
	 * @return The swaprates of the Bermudan swaption.
	 */
	public double[] getSwaprates() {
		return swaprates;

	}
}
