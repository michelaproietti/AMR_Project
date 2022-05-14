#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>

namespace ifopt {
using Eigen::Vector2d;

Eigen::Vector3d qIni, qFin;

// tunable
double w_v = 0.5;
double w_omega = 0.5;
double w_delta = 0.5;
int N = 5;
double vMin = -1.0;
double vMax = 1.0;
double omegaMin = -2.0;
double omegaMax = 2.0;
double deltaMin = 0.0;
double deltaMax = 100.0;

// not tunable
int numDecVar = 3*(N+1)+3*N;
int numConstr = 3*N+6;
	
void setTPBVP(Eigen::Vector3d _qIni, Eigen::Vector3d _qFin){
	qIni = _qIni;
	qFin = _qFin;	
}

class ExVariables : public VariableSet {
	public:
  
	ExVariables() : ExVariables("var_set1") {};
	ExVariables(const std::string& name) : VariableSet(numDecVar, name) 
	{
		// the initial values where the NLP starts iterating from
		decVar.resize(numDecVar);	
		decVar.setZero();
	}

	void SetVariables(const VectorXd& x) override  
	{
		decVar = x;
	};

	VectorXd GetValues() const override  
	{
		VectorXd val(numDecVar);
		val = decVar;
		return val;
	};

	// Each variable has an upper and lower bound set here
	VecBound GetBounds() const override 
	{
		VecBound bounds(GetRows());
		for(int i = 0; i < 3*(N+1); i++) bounds.at(i) = NoBound;
		for(int i = 3*(N+1); i < 3*(N+1)+N; i++) bounds.at(i) = Bounds(vMin, vMax); 
		for(int i = 3*(N+1)+N; i < 3*(N+1)+N+N; i++) bounds.at(i) = Bounds(omegaMin, omegaMax); 
		for(int i = 3*(N+1)+N+N; i < 3*(N+1)+N+N+N; i++) bounds.at(i) = Bounds(deltaMin, deltaMax); 

		return bounds;
	}

	private:
		VectorXd decVar;
};


class ExConstraint : public ConstraintSet {
	public:
	ExConstraint() : ExConstraint("constraint1") {}

	ExConstraint(const std::string& name) : ConstraintSet(numConstr, name) {}

	VectorXd GetValues() const override 
	{
		VectorXd g(GetRows());
		VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();

		for(int k = 0; k < N; k++){
			double x_k = x(k);
			double x_kp1 = x(k+1);
			double y_k = x(k+N+1);
			double y_kp1 = x(k+N+2);
			double theta_k = x(k+2*(N+1));
			double theta_kp1 = x(k+2*(N+1)+1);
			double v_k = x(k+3*(N+1));
			double omega_k = x(k+3*(N+1)+N);
			double delta_k = x(k+3*(N+1)+2*N);
			
			g(3*k) = x_kp1 - x_k - delta_k * v_k * cos(theta_k);
			g(3*k+1) = y_kp1 - y_k - delta_k * v_k * sin(theta_k);
			g(3*k+2) = theta_kp1 - theta_k - delta_k * omega_k;		
		} 

		double x_0 = x(0);
		double y_0 = x(N+1);
		double theta_0 = x(2*N+2);
		double x_N = x(N);
		double y_N = x(2*N+1);
		double theta_N = x(3*N+2);
			
		g(3*N) = x_0;
		g(3*N+1) = y_0;
		g(3*N+2) = theta_0;
		g(3*N+3) = x_N;
		g(3*N+4) = y_N;
		g(3*N+5) = theta_N;

		return g;
	};

	VecBound GetBounds() const override  
	{
		VecBound b(GetRows());
    
		for(int k = 0; k < N; k++){
			b.at(3*k) = Bounds(0.0, 0.0);
			b.at(3*k+1) = Bounds(0.0, 0.0);
			b.at(3*k+2) = Bounds(0.0, 0.0);
		} 

		b.at(3*N) = Bounds(qIni(0), qIni(0));
		b.at(3*N+1) = Bounds(qIni(1), qIni(1));
		b.at(3*N+2) = Bounds(qIni(2), qIni(2));
		b.at(3*N+3) = Bounds(qFin(0), qFin(0));
		b.at(3*N+4) = Bounds(qFin(1), qFin(1));
		b.at(3*N+5) = Bounds(qFin(2), qFin(2));

		return b;
	}

	void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override
	{
		if (var_set == "var_set1") {
			// we use the approximate Jacobians so we don't override this function
		}
	}
};


class ExCost: public CostTerm {
	public:
	ExCost() : ExCost("cost_term1") {}
	ExCost(const std::string& name) : CostTerm(name) {}

	double GetCost() const override 
	{
		VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();

		double cost_v = 0.0;
		double cost_omega = 0.0;
		double cost_delta = 0.0;

		for(int i = 3*(N+1); i < 3*(N+1)+N; i++) cost_v = x(i)*x(i); 
		for(int i = 3*(N+1)+N; i < 3*(N+1)+N+N; i++) cost_omega = x(i)*x(i); 
		for(int i = 3*(N+1)+N+N; i < 3*(N+1)+N+N+N; i++) cost_delta = x(i)*x(i); 

		double cost = 1 + w_v * cost_v + w_omega * cost_omega + w_delta * cost_delta;

		return cost;
	};

	void FillJacobianBlock (std::string var_set, Jacobian& jac) const override 
	{
	    if (var_set == "var_set1") {
			VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();

			for(int i = 0; i < 3*(N+1); i++) jac.coeffRef(0, i) = 0.0; // derivative of cost w.r.t qi 
			for(int i = 3*(N+1); i < 3*(N+1)+N; i++) jac.coeffRef(0, i) = 2.0 * w_v * x(i);	// derivative of cost w.r.t vi 
			for(int i = 3*(N+1)+N; i < 3*(N+1)+N+N; i++) jac.coeffRef(0, i) = 2.0 * w_omega * x(i);	// derivative of cost w.r.t omegai
			for(int i = 3*(N+1)+N+N; i < 3*(N+1)+N+N+N; i++) jac.coeffRef(0, i) = 2.0 * w_delta * x(i); // derivative of cost w.r.t deltai
		}
	}
};

} // namespace opt




































