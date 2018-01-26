struct Observs {
	vector<Est> obs;
	Sensor sen;

	Observs () {
	};

	Observs (Sensor sen) : sen(sen) {
	}

	Observs (vector<Est> obs, Sensor sen) : obs(obs), sen(sen) {
	}
};