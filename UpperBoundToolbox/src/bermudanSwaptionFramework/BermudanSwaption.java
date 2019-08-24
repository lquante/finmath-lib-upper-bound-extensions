package bermudanSwaptionFramework;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;

/**
 * Implements a Bermudan swaption product..
 * 
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
	 * Method using the valuation method specified in the option.
	 */
	@Override
	public RandomVariable getValue(double evaluationTime, LIBORModelMonteCarloSimulationModel model)
			throws CalculationException {
		return valuationMethod.getValueEstimation(this, evaluationTime, model, null);
	}

	/**
	 * Method to get a clone with another valuation method
	 * 
	 * @param valuationMethod The valuation method to be used.
	 * @return A Bermudan Swaption with the input valuation method
	 */
	public BermudanSwaption getBermudanSwaptionWithChangedValuationMethod(
			BermudanSwaptionValueEstimatorInterface valuationMethod) {
		return new BermudanSwaption(this.getCurrency(), this.getIsPeriodStartDateExerciseDate(), this.getFixingDates(),
				this.getPeriodLengths(), this.getPaymentDates(), this.getPeriodNotionals(), this.isCallable(),
				this.swaprates, valuationMethod);
	}

	/**
	 * Shifts the first exercise date of an Bermudan swaption needed e.g. in method
	 * of Andersen-Broadie.
	 * 
	 * @param startingIndex first fixing date to be kept
	 * @return Clone of the BermudanSwaption with the input starting index.
	 */
	public BermudanSwaption getCloneWithModifiedStartingPeriod(int startingIndex) {

		if (startingIndex >= this.getFixingDates().length)
			startingIndex = this.getFixingDates().length - 1;

		// shorten all arrays for a new BermudanSwaption

		return getCloneWithModifiedStartingAndFinalPeriod(startingIndex, this.getFixingDates().length - 1);

	}

	/**
	 * Shifts the first and final exercise date of an Bermudan swaption needed e.g.
	 * in method of Andersen-Broadie.
	 * 
	 * @param startingIndex first fixing date to be kept
	 * @param finalIndex    last fixing date to be kept (inclusive)
	 * @return Clone of the BermudanSwaption with the input starting index and final
	 *         incex.
	 */
	public BermudanSwaption getCloneWithModifiedStartingAndFinalPeriod(int startingIndex, int finalIndex) {

		if (startingIndex >= this.getFixingDates().length)
			startingIndex = this.getFixingDates().length - 1;
		finalIndex += 1;
		if (finalIndex >= this.getFixingDates().length)
			finalIndex = this.getFixingDates().length;

		// shorten all arrays for a new BermudanSwaption

		double[] adjustedFixingDates = copyDoubleArray(this.getFixingDates(), startingIndex, finalIndex);
		boolean[] adjustedIsPeriodStartDateExerciseDate = copyBooleanArray(this.getIsPeriodStartDateExerciseDate(),
				startingIndex, finalIndex);

		double[] adjustedPaymentDates = copyDoubleArray(this.getPaymentDates(), startingIndex, finalIndex);
		double[] adjustedPeriodLengths = copyDoubleArray(this.getPeriodLengths(), startingIndex, finalIndex);
		double[] adjustedPeriodNotionals = copyDoubleArray(this.getPeriodNotionals(), startingIndex, finalIndex);
		double[] adjustedSwapRates = copyDoubleArray(this.getSwaprates(), startingIndex, finalIndex);

		return new BermudanSwaption(this.getCurrency(), adjustedIsPeriodStartDateExerciseDate, adjustedFixingDates,
				adjustedPeriodLengths, adjustedPaymentDates, adjustedPeriodNotionals, this.isCallable(),
				adjustedSwapRates, this.getValuationMethod());

	}

	// methods to copy a part of an array using System.arraycopy

	double[] copyDoubleArray(double[] arrayToBeCopied, int startingIndex, int finalIndex) {
		double[] copiedArray = new double[finalIndex - startingIndex];
		System.arraycopy(arrayToBeCopied, startingIndex, copiedArray, 0, finalIndex - startingIndex);
		return copiedArray;

	}

	boolean[] copyBooleanArray(boolean[] arrayToBeCopied, int startingIndex, int finalIndex) {
		boolean[] copiedArray = new boolean[finalIndex - startingIndex];
		System.arraycopy(arrayToBeCopied, startingIndex, copiedArray, 0, finalIndex - startingIndex);
		return copiedArray;

	}

	// some getters

	/**
	 * @return the valuation method of the Bermudan swaption
	 */
	public BermudanSwaptionValueEstimatorInterface getValuationMethod() {
		return valuationMethod;
	}

	/**
	 * Changes the valuation method of an Bermudan swaption without generating a new
	 * object.
	 * 
	 * @param valuationMethod The valuation method to be used.
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
