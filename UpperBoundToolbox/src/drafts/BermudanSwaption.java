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
	 * @param startingDate timepoint at which the option should restart (first
	 *                     implementation: only from the exisiting fixing dates)
	 * @return
	 */
	public BermudanSwaption getCloneWithModifiedStartingPeriod(double startingDate) {
		int startingIndex = 0;
		// determine first fixing date to be used
		for (int i = 0; i < this.getFixingDates().length; i++) {
			if (this.getFixingDates()[i] == startingDate) {
				startingIndex = i;
				break;
			}
		}

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

	public BermudanSwaption getCloneWithModifiedFinalPeriod(double finalDate) {
		int finalIndex = this.getPaymentDates().length - 1;
		// determine last payment date to be used
		for (int i = this.getPaymentDates().length - 1; i >= 0; i--) {
			if (this.getPaymentDates()[i] == finalDate) {
				finalIndex = i + 1;
				break;
			}
		}

		// shorten all arrays for a new BermudanSwaption

		double[] adjustedFixingDates = Arrays.copyOfRange(this.getFixingDates(), 0, finalIndex);
		boolean[] adjustedIsPeriodStartDateExerciseDate = Arrays.copyOfRange(this.getIsPeriodStartDateExerciseDate(), 0,
				finalIndex);
		double[] adjustedPeriodLengths = Arrays.copyOfRange(this.getPeriodLengths(), 0, finalIndex);
		double[] adjustedPaymentDates = Arrays.copyOfRange(this.getPaymentDates(), 0, finalIndex);
		double[] adjustedPeriodNotionals = Arrays.copyOfRange(this.getPeriodNotionals(), 0, finalIndex);
		double[] adjustedSwapRates = Arrays.copyOfRange(this.getSwaprates(), 0, finalIndex);
		BermudanSwaptionValueEstimatorInterface unadjustedValuationMethod = this.getValuationMethod();
		return new BermudanSwaption(this.getCurrency(), adjustedIsPeriodStartDateExerciseDate, adjustedFixingDates,
				adjustedPeriodLengths, adjustedPaymentDates, adjustedPeriodNotionals, this.isCallable(),
				adjustedSwapRates, unadjustedValuationMethod);

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
