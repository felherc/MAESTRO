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

import java.util.ArrayList;

import maestro.solution.Solution;

public class RosenbrockSol implements Solution
{
	
	private String id;
	
	private double x;
	
	private double y;
	
	private double f;

	public RosenbrockSol(String id, double x, double y)
	{		
		this.id = id;
		this.x = x;
		this.y = y;

		f = (1 - x)*(1 - x) + 100*(y - x*x)*(y - x*x);
		//f = Math.sqrt((750 - x)*(750 - x) + (750 - y)*(750 - y));
	}
	
	public double getX()
	{
		return x;
	}
	
	public double getY()
	{
		return y;
	}
	
	public double getF()
	{
		return f;
	}
	
	@Override
	public int compareTo(Solution otherSol) 
	{
		RosenbrockSol other = (RosenbrockSol) otherSol;
		if(x == other.getX() && y == other.getY())
			return 0;
		else if(f < other.getF())
			return 1;
		else
			return -1;
	}

	@Override
	public Solution createNew(int id, ArrayList<Integer> discValues,
			ArrayList<Double> contValues) 
	{
		String ids = "Solution " + id;
		RosenbrockSol sol = new RosenbrockSol(ids, contValues.get(0), contValues.get(1));
		return sol;
	}

	@Override
	public String getId() 
	{
		return id;
	}

	@Override
	public ArrayList<Integer> getDiscValues() 
	{
		return null;
	}

	@Override
	public ArrayList<Double> getContValues() 
	{
		ArrayList<Double> values = new ArrayList<Double>();
		values.add(x);
		values.add(y);
		return values;
	}

	@Override
	public String getReportHeader() 
	{
		return "f(x, y)";
	}

	@Override
	public String getReport() 
	{
		return "" + f;
	}

	@Override
	public boolean optimizationConverged() 
	{
		return false;
	}
	
}
