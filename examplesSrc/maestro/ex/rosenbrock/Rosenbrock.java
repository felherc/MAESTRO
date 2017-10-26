/**
 * Copyright 2017 Felipe Hernández
 * 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License. You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software distributed under the 
 * License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 * express or implied. See the License for the specific language governing permissions and
 * limitations under the License.
 */

package maestro.ex.rosenbrock;

import java.io.IOException;

import maestro.MAESTROptimizer;
import maestro.Monitor;

public class Rosenbrock implements Monitor 
{

	private MAESTROptimizer maestro;
	
	public static void main(String[] args) 
	{
		Rosenbrock rosenbrock = new Rosenbrock();
		rosenbrock.run();
	}
	
	public void run()
	{
		RosenbrockSol solution = new RosenbrockSol("Default solution", 0, 0);
		maestro = new MAESTROptimizer(solution, this, false, 20);
		
		maestro.addContVar("x", -10, 10);
		maestro.addContVar("y", -10, 10);
		maestro.startOptimization(2*60*1000, 100000);
	}

	@Override
	public void terminate(String message) 
	{
		try 
		{
			maestro.writeReport("data/tests/Rosenbrock results.txt", true, true, false);
		} catch (IOException e) 
		{
			e.printStackTrace();
		}
	}

	@Override
	public void reset() 
	{
		// Do nothing
	}

}
