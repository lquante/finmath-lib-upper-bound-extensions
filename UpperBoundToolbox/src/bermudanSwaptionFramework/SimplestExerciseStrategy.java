package bermudanSwaptionFramework;

import net.finmath.stochastic.RandomVariable;

public class SimplestExerciseStrategy implements ExerciseStrategyInterface {

	/** 
	 * Implements the simplest exercise strategy: Exercise if the exercise value is bigger or equal the continuation value.
	 * @see drafts.ExerciseStrategyInterface#getExerciseIndicators(net.finmath.stochastic.RandomVariable, net.finmath.stochastic.RandomVariable)
	 */
	@Override
	public RandomVariable getExerciseIndicators(RandomVariable continuationValue, RandomVariable exerciseValue) {

		return exerciseValue.sub(continuationValue);
	}

}
