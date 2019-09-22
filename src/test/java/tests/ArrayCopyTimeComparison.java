package tests;

import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;

/**	This class tests if the chosen array copy method using the System library is indeed faster than the alternative via Array library.
 * @author Lennart Quante
 *  @version 1.0
 */	

public class ArrayCopyTimeComparison {
	
	

	@Test
	public void compareArrayCopyTimes() {

		int numberOfRepetitions = 100;
		int startingIndex = 1;
		int finalIndex = 8;

		double[] testArray = { 1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9, 9.10, 10.11 };

		for (int magnitudeIndex = 5; magnitudeIndex < 15; magnitudeIndex++) {
			numberOfRepetitions = (int) Math.pow(10, magnitudeIndex);

			long timingArrayCopyOfRangeStart = System.currentTimeMillis();

			for (int i = 0; i < numberOfRepetitions; i++) {
				double[] copiedArray = Arrays.copyOfRange(testArray, startingIndex, finalIndex);
			}
			long timingArrayCopyOfRangeEnd = System.currentTimeMillis();

			long timeArrayCopyOfRange = timingArrayCopyOfRangeEnd - timingArrayCopyOfRangeStart;

			long timingSystemCopyStart = System.currentTimeMillis();

			for (int i = 0; i < numberOfRepetitions; i++) {
				double[] copiedArray = new double[finalIndex - startingIndex];
				System.arraycopy(testArray, startingIndex, copiedArray, 0, finalIndex - startingIndex);
			}
			long timingSystemCopyEnd = System.currentTimeMillis();
			long timeSystemCopy = timingSystemCopyEnd - timingSystemCopyStart;
			System.out.println("Time comparison for a repeated array copying, repetitions:" + numberOfRepetitions);
			System.out.println("Time needed using arrays.copyOfRange: " + timeArrayCopyOfRange);
			System.out.println("Time needed using SystemCopy: " + timeSystemCopy);
			Assert.assertTrue(timeArrayCopyOfRange>timeSystemCopy);
		}
	}
}
