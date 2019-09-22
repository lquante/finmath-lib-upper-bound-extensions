package simulationMethods;

import bermudanswaptionframework.BermudanSwaption;
import bermudanswaptionframework.BermudanSwaptionValueEstimatorInterface;

/**
 * Class to construct a factory for test swaptions.
 * @author Lennart Quante
 *  @version 1.0
 */
public class TestSwaptionFactory {

	private int numberOfExercisePeriods;

	/**
	 * Constructs a factory to provide swaptions with the following parameters fixed:
	 * 
	 * @param numberOfExercisePeriods number of exercise periods
	 * @param periodLength length of option periods
	 * @param swaprate swaprate to be used
	 */
	public TestSwaptionFactory(int numberOfExercisePeriods, double periodLength, double swaprate) {
		super();
		this.numberOfExercisePeriods = numberOfExercisePeriods;
		this.periodLength = periodLength;
		this.swaprate = swaprate;
	}

	private double periodLength;
	private double swaprate;

	/**
	 * Method to construct callable Bermudan swaption for testing purposes from
	 * reduced input. All fixing dates are exercise dates, equal period length,
	 * notional ==1 and constant swaprate
	 * 
	 * @param firstFixingDate         first fixing date of the swaption
	 * @param valuationMethod         valuation method to be used
	 * @return A Bermudan Swaption with the given parameters
	 */
	public BermudanSwaption constructBermudanSwaption(double firstFixingDate,
			BermudanSwaptionValueEstimatorInterface valuationMethod) {

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
		return new BermudanSwaption("EURO", isPeriodStartDateExerciseDate, fixingDates, periodLengths, paymentDates,
				notionals, true, swaprates, valuationMethod);

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
