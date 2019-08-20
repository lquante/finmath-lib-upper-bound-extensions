package bermudanSwaptionFramework;

import net.finmath.stochastic.RandomVariable;

/**
  Implementing the basic exercise strategy: Exercise if the exercise value is bigger or equal the continuation value.
 * @author Lennart Quante
 * @version 1.0
 */
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