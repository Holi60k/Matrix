#include <iostream>

template<class T>
class Matrix {

public:

	
	Matrix(int a = 0, int b = 0):siX(a),siY(b)
	{
		//inicializálás itt történik, memóriában foglalunk helyet ennek a sok szép számnak...
		//beállítjuk a mutatónkat hogy mutasson egy mutatóra... dafaq... igen ez most így lesz :)
		Matx = new T*[a];
		for(int i = 0; i < a; i++) {
			//aztán pedig egyesével lefoglalunk a mutatók számára dinamikus memóriát
			Matx[i] = new T[b];
		}

		if(siX == siY) {

			Squared = true;
		} else {
			Squared = false;
		}
	};

	~Matrix()
	{
		//destruktor törli a sok szép számunkat
		//logikus hogy előbb azokat az adatokat töröljük amiket legutoljára raktunk a memóriába
		//pointer to pointer, A to B... B(k) törlése előbb a memóriából aztán pedig az A-t
		for(int i = 0; i < this->GetX(); i++)
		{
			delete [] Matx[i];
		
		}
		delete [] Matx;

	};
	//értékadó operátorunk, másoló értékadás!
	Matrix & operator= (Matrix & A) {

		//std::cout << "copy assingment" << std::endl;
		siX = A.siX;
		siY = A.siY;
		//egy előre lefoglalt kibaszott 0 terület törlése :'(
		delete [] Matx;


		Matx = new T*[siX];
		for(int i = 0; i < siX; i++) {
			//aztán pedig egyesével lefoglalunk a mutatók számára dinamikus memóriát
			Matx[i] = new T[siY];
		}

		for(int i = 0; i < A.GetX();i++) {
			for(int j = 0; j < A.GetY();j++) {
				this->FillMatrix(A.GetValue(i,j),i,j);
			}
		}

		Squared = A.Squared;
		return *this;
			
	}

	Matrix (Matrix & t) {

		//std::cout << "copy ctor " << std::endl;
		siX = t.GetX();
		siY = t.GetY();

		Matx = new T*[siX];
		//std::cout << "Old: " << &t.Matx << std::endl;
		//std::cout << "New: " << &Matx << std::endl;
		for(int i = 0; i < siX; i++) {
			//aztán pedig egyesével lefoglalunk a mutatók számára dinamikus memóriát
			Matx[i] = new T[siY];
			//std::cout << "Old: "<< i << " - " << &t.Matx[i] << std::endl;
			//std::cout << "New: "<< i << " - " << &Matx[i] << std::endl;
		}

		for(int i = 0; i < t.GetX();i++)
			for(int j = 0; j < t.GetY();j++) {
				FillMatrix(t.GetValue(i,j),i,j);
			}

		//GetWholeMatrix();
		Squared = t.Squared;

	}
	//Matrix (const Matrix&);
	//két mátrix összeadására szolgáló operátorunk
	Matrix & operator+ (const Matrix & t) {

		std::cout << "Starting add the two Matrixes..." << std::endl;
		if(this->GetY() == t.GetX()) {

			Matrix<T> Result(this->GetX(),t.GetY());

			for(int i = 0; i < Result.GetX();i++)
			{
				for(int j = 0; j < Result.GetY(); j++)
				{				
					Result.FillMatrix(this->GetValue(i,j)+t.GetValue(i,j),i,j);
				}

			}

			*this = Result;
			return *this;

		} else {
				std::cout << "Sorry I can not multiplicate these two Matrixes..."<< std::endl;
				
		}
		 return *this;
	}

	Matrix & operator* (const Matrix & A)
	{
		T Var = 0;
		std::cout << "Starting multiplicate the two Matrixes..." << std::endl;
		if(this->GetY() == A.GetX())
		{
			Matrix<T> Result(this->GetX(),A.GetY());
			
			for(int i = 0; i < Result.GetX();i++)
			{
				for(int j = 0; j < Result.GetY(); j++)
				{				
					
					for(int l = 0; l < this->GetY();l++ )
					{
						Var += this->GetValue(j,l) * A.GetValue(l,i);
						//std::cout << "Actual Value:" << Var << "\nthis("<<j<<","<<l<<"):"<<this->GetValue(j,l)<<"\nA("<<l<<","<<i<<"):"<<A.GetValue(l,i)<<std::endl;					
					}
					//std::cout << "Actual Value:" << Var << "\nthis:"<<this->GetValue(j,i)<<"\nA:"<<A.GetValue(j,i)<<std::endl;					
					Result.FillMatrix(Var,j,i);
					Var = 0;
					//std::cout << "Actual Value:" << Var <<"\n"<<std::endl;					

				}

			}

			*this = Result;
			return *this;
		}
		else
		{
			std::cout << "Sorry I can not multiplicate these two Matrixes..."<< std::endl;
			
		}
	}
	//Mátrix feltöltése, v - érték, x - x koordináta, y - y koordináta
	void FillMatrix(T v, int x, int y) const {
		Matx[x][y] = v;		
	}

	int GetX() const {
		return siX;
	}

	int GetY() const {
		return siY;
	}

	T GetValue(int x, int y) const {
		return Matx[x][y];
	}

	void ChangeRows(int a, int b) {
		T *cs;
		cs = Matx[a];
		Matx[a] = Matx[b];
		Matx[b] = cs;
		Det_Sign *= -1;
		
		
	}

	T* GetRows(int i) { 
		//std::cout<< i << " - " << Matx[i] << std::endl;
		return Matx[i];
	}
	int Get_Det_Sign() {
		return Det_Sign;
	}
	void GaussElimination() {
		std::cout << "Gauss" << std::endl;

		for(int i = 0; i < siX; i++) {
			ChangeRows(0,i);
			//std::cout << "Multiplicate - " << i << "\nOsztok:" << 1/Matx[0][i]<<std::endl;

			MultiplicateRow(0,Matx[0][i] != 0?1/Matx[0][i]:1);
			//this->GetWholeMatrix();
			//std::cout << "-" << std::endl;

			//std::cout << "AddRow2Row - " << i+1 <<std::endl;
			for(int j = 1; j < siX; j++) {

				AddRow2Row(j,0,(-1)*Matx[j][i]);
			//	this->GetWholeMatrix();
			//	std::cout << "--" << std::endl;
			}
			ChangeRows(0,i);
			//std::cout << "---" << std::endl;
			//this->GetWholeMatrix();
			//std::cout << "--------------" << std::endl;

			
			
		}

	}
	//ez inkább felső 3-szög alakra hozza a mátrixot :)
	void GaussElimination_2() {
		//std::cout << "Gauss_2" << std::endl;
		T pivot;
		for(int i = 0; i < siX; i++) {
			if(WhereIsItsRow(i) == i){
				pivot = Matx[i][i] != 0?1/Matx[i][i]:1;
				MultiplicateRow(i,pivot);

				for(int j = i+1; j < siX; j++) {

					AddRow2Row(j,i,(-1)*Matx[j][i]);
				}

				MultiplicateRow(i,1/pivot);
			}
			else {
				ChangeRows(i,WhereIsItsRow(i));
				GaussElimination_2();
			}
	}

		//this->GetWholeMatrix();

	}
	void AddRow2Row(int a, int b, T value) {

		for(int i = 0; i < siY; i++) {

			Matx[a][i] += value * Matx[b][i];
		}

	}
	void GaussEliminationWVector(Matrix & b) {

		std::cout << "Gauss_W_Vector" << std::endl;
		T pivot;
		int Free_Param = 0;
		for(int i = 0; i < siX; i++) {
			ChangeRows(0,i);
			b.ChangeRows(0,i);
			pivot = Matx[0][i] != 0?1/Matx[0][i]:1;
			
			MultiplicateRow(0,pivot);

			b.MultiplicateRow(0,pivot);
			//this->GetWholeMatrixVector(b);
			for(int j = 1; j < siX; j++) {
				
				pivot = (-1)* Matx[j][i];
				
				b.AddRow2Row(j,0,pivot);
				AddRow2Row(j,0,pivot);
			}

		
			ChangeRows(0,i);
			b.ChangeRows(0,i);
	
			
		}
		Free_Param = siY - Rank();
		std::cout << "Free Parameters: " << Free_Param<< std::endl;

	}

	int Determinant() {

		//std::cout << "Determinant" << std::endl;
		int Det = 1;
		Matrix<T> Result = *this;
		if(Squared) {

			if(!Result.IsTriangle()) {
				if(Result.CanIReOrderToTriangle()) 
					Result.ReOrderTriangle();
			}
			
			Result.GaussElimination_2();
			for(int i = 0; i < siX; i++) {
					Det *= Result.Matx[i][i];
			}
		
			return Det*Result.Get_Det_Sign();

		} else {
			std::cout << "Sorry I can not count its determinant cause this matrix is non-squared." << std::endl;
			return -1;
		}
	}

	int Rank() {
		Matrix<T> Result = *this;
		int rank = 0;
		Result.GaussElimination();

		for(int i = 0; i < siX; i++) {

			if(Result.CountZeroItems(i) == 0) rank++;
			std::cout << "R: " << rank << std::endl;
		}

		return rank;

	}
	int CountZeroItems(int a) {
		int Zero = 0;

		if(Matx[a][a] != 0 && Squared) {
			return 0;
		}

		for(int i = 0; i < siY; i++) {
			if(Matx[a][i] == 0) Zero++;
		}

		//std::cout<< "Z: " << a << " - " << Zero << std::endl;
		//std::cout <<"A: " << a << " - " << Matx[a][a] << std::endl;
		
		return Zero;
	}
	void ReOrderRows() { 
		int Zero_Num[siX] = {0};

		for(int i = 0; i < siX; i++) {
		Zero_Num[i] = 0;
			for(int j= 0; j < siY; j++) {
				if(Matx[i][j] == 0)
					Zero_Num[i]++;
			}
		}
	}

	bool IsFineRow(int a) {

		int Zero_Num = 0;

		for(int i = 0; i < a; i++) {
			if(Matx[a][i] == 0)
				Zero_Num++;
		}
		if(Zero_Num == a)
			return true;
		else
			return false;

	}

	int WhereIsItsRow(int a) {

		int Zero_Num = 0;
			for(int j = 0; j < a; j++) {
				
				if(Matx[a][j] == 0) {
					Zero_Num++;
				}
			}
		return Zero_Num;

	}
	void ReOrderTriangle(){
		// 0 0 1
		// 1 0 1
		// 0 1 1

		// 1 0 1
		// 0 1 1
		// 0 0 1

		int Places[siX];
		int Changes = 0;
		for(int i = 0; i < siX; i++) {
			Places[i] = i;
			if(!IsFineRow(i)) {
				ChangeRows(i,WhereIsItsRow(i));
				Places[i] = WhereIsItsRow(i);
			} 
			
		}
	}

	bool IsTriangle() {

		int Zero_Num[siX];
		bool Row_OK[siX];
		int Rows = 0;

		for(int i = 0; i < siX; i++) {
		Row_OK[i] = true;
			for(int j = 0; j < i; j++) {
				if(Matx[i][j] == 0) {
					Row_OK[i] = true;
				} else {
					Row_OK[i] = false;
					break;
				}
			}
		}

		for(int i = 0; i < siX; i++){
			if(Row_OK[i] == true) {
				Rows++;
				
			}
		}

		if(Rows == siX)
			return true;
		else
			return false;

	}

	bool CanIReOrderToTriangle() {
		int Zero_Num[siX];
		
		int Rows = 0;
		int a = 0;
		for(int i = 0; i < siX; i++) {
		Zero_Num[i] = 0;
			while(Matx[i][a] == 0 && a < siY) {
				//std::cout << Matx[i][a] << " ";
				Zero_Num[i]++;
				a++;
			}
			//std::cout << std::endl;
			a = 0;
			//std::cout <<"Folytonos nullák az "<< i << " sorban:" << Zero_Num[i] << std::endl;
		}
		for(int i = 0; i < siX; i++) {
		
			if(Zero_Num[i] > 0){
				Rows++;
			}

		}
		if(Rows == siX-1)
			return true;
		else
		 return false;

	}

	void MultiplicateRow(int a, T Multi) {

		for(int i = 0; i < siY; i++) {
			//std::cout << "SubRow:(in) " << i  << " - " << Matx[a][i] << std::endl;
			Matx[a][i] *= Multi;

		}
	}

	void GetWholeMatrix() const
	{
		std::cout << "---" << std::endl;
		//std::cout << "Get the whole Martix..." << std::endl;
		for(int i = 0; i < siX; i++)
		{
			for(int j = 0; j < siY; j++)
			{
				//Matx[i][j] = v;
				std::cout << Matx[i][j] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "---" << std::endl;

	}
	void GetWholeMatrixVector(Matrix & b) const
	{
		//std::cout << "Get the whole Martix..." << std::endl;
		for(int i = 0; i < siX; i++)
		{
			for(int j = 0; j < siY; j++)
			{
				//Matx[i][j] = v;
				std::cout << Matx[i][j] << " ";
			}
			std::cout << b.Matx[i][0];
			std::cout << std::endl;
		}

	}
	Matrix Zero_Vector() {

	}
private:
	int siX;
	int siY;
	bool Squared;
	bool Null_Matrix;
	int Det_Sign = 1;
	T** Matx = NULL;
};



int main()
{
	float x = 3,y = 3;
	//std::cout << "How many rows do you want:";
	//std:: cin >> x;
	//std::cout << "How many collumns do you want:";
	//std::cin >> y;
	
	Matrix<float> A(x,y), W(x,1);

	int v=0;
	std::cout << "Fill the first Matrix" << std::endl;
	/*for(int i = 0; i < x; i++)
	{
		for(int j = 0; j<y; j++)
		{
			std::cin >> v;
			//v++;
			A.FillMatrix(v,i,j);
		}
	}*/
	//X-ek
	A.FillMatrix(0,0,0);
	A.FillMatrix(1,0,1);
	A.FillMatrix(1,0,2);

	//Y-ok
	A.FillMatrix(1,1,0);
	A.FillMatrix(1,1,1);
	A.FillMatrix(0,1,2);
	//Z-k
	A.FillMatrix(0,2,0);
	A.FillMatrix(0,2,1);
	A.FillMatrix(1,2,2);
	
	std::cout << "Fill the Vector" << std::endl;
	W.FillMatrix(2,0,0);
	W.FillMatrix(-2,1,0);
	W.FillMatrix(-3,2,0);
	/*for(int i = 0; i < x; i++)
	{
		for(int j = 0; j<1; j++)
		{
			std::cin >> v;
			//v++;
			W.FillMatrix(v,i,j);
		}
	}*/
	//std::cout << "Fill the second Matrix" << std::endl;
	/*for(int i = 0; i < x; i++)
	{
		for(int j = 0; j<y; j++)
		{
			std::cin >> v;
			B.FillMatrix(v,i,j);
		}
	}*/
	A.GetWholeMatrix();
	//B.GetWholeMatrix();
	//A.ChangeRows(0,1);
	//A.GetWholeMatrix();
	//A.GaussEliminationWVector(W);
    //A.GetWholeMatrixVector(W);
	//std::cout << "A Determinant:" << A.Determinant() << std::endl;
	if(A.IsTriangle()) {
		std::cout << "YEAH it is a triangle formed!" << std::endl;
		//std::cout << "A1 Determinant:" << A.Determinant() << std::endl;
	} else {
		if(A.CanIReOrderToTriangle())
			//A.ReOrderTriangle();
		std::cout << "A2 Determinant:" << A.Determinant() << std::endl;
	}
	
	//A.GetWholeMatrix();
	//A.GaussElimination_2();
	A.GetWholeMatrix();
	//std::cout << "A rangja:" << A.Rank() << std::endl;
	//B.GetWholeMatrix();
	//C = A;
	//C.GetWholeMatrix();
	//amikor meghívódik az operator()* fgv akkor this = A...
	//C = A+B;
	//C.GetWholeMatrix();
	return 0;
}