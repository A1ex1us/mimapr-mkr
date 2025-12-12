#include <algorithm>
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <fstream>

class Form { //Класс фигуры
protected:
    static size_t counter_;
    size_t id_;
    bool excluded_;
    int _boundtype;

public:
    Form();

    virtual double Function(double, double);

    virtual std::pair<double, double> Deriative(double, double);

    virtual bool Inhere(double, double);
    
    virtual std::pair<double, double> missX(double); 
    virtual std::pair<double, double> missY(double); 
    virtual std::pair<double, double> size();
    virtual int GetB();
    [[nodiscard]] size_t Get_ID() const;

    [[nodiscard]] bool Excluded() const;

    bool operator==(size_t) const;
};

class Rectangle : public Form { //Класс прямоугольника
private:
    double a_;
    double b_;
    double h_x_;
    double h_y_;
public:
    Rectangle(double, double, double, double, bool, int);

    double Function(double, double) override;

    std::pair<double, double> Deriative(double, double) override;

    bool Inhere(double, double) override;
    std::pair<double, double> missX(double) override;
    std::pair<double, double> missY(double) override;
    std::pair<double, double> size() override;
};

class Circle : public Form {  //Класс окружности
private:
    double a_;
    double b_;
    double h_x_;
    double h_y_;
public:
    Circle(double, double, double, double, bool, int);

    double Function(double, double) override;

    std::pair<double, double> Deriative(double, double) override;

    bool Inhere(double, double) override;
    std::pair<double, double> missX(double) override;
    std::pair<double, double> missY(double) override;
    std::pair<double, double> size() override;
};

class Arc : public Form {  //Класс дуги
private:
    double a_;
    double b_;
    double h_x_;
    double h_y_;
public:
    Arc(double, double, double, double, bool, int);

    double Function(double, double) override;

    std::pair<double, double> Deriative(double, double) override;
    std::pair<double, double> missX(double) override;
    std::pair<double, double> missY(double) override;
    std::pair<double, double> size() override;
    bool Inhere(double, double) override;
};

class Object {   //Класс объектов
private:
    std::vector<Form*> forms_;
    double _w;
    double _h;
    void Updsize();
public:
    Object();
    ~Object();

    double Inhere(double, double);
    double Width() const;
    double Height() const;
    bool Add_Form(const std::string&, std::map<std::string, double>&, bool, int);
   

    bool Delete_Form(size_t);
    std::pair<double, double> Filly(double, double);
    std::pair<double, double> Fillx(double, double);

    std::vector<size_t> Get_IDs();
    Form* Who(double, double);
};

class Node   //Класс элементов сетки
{
    double _x;
    double _y;
    double _t;
    int _btype;
    
    Node* _left;
    Node* _right;
    Node* _above;
    Node* _below;

    double _node_H_robin;       // Персональный H = h/k для этого узла
    double _node_T_ambient_robin; 
public:
    Node(double x = 0., double y = 0., int type = 0., double t = 0.)
        : _x(x), _y(y), _t(t), _left(nullptr), _right(nullptr), _above(nullptr), _below(nullptr),
          _btype(type), _node_H_robin(0.0), _node_T_ambient_robin(0.0) {} // Инициализация

    double T() const;
    double X() const;
    double Y() const;
    double Dist(const Node*) const;
    void LinkX(Node*, Node*);
    void LinkY(Node*, Node*);
    Node*& l();
    Node*& r();
    Node*& u();
    Node*& d();
    void SetT(double);
    bool IsBound();
    void SetB(int);

    void SetRobinParameters(double H, double T_ambient) {
        _node_H_robin = H;
        _node_T_ambient_robin = T_ambient;
    }
    // Вспомогательный метод для получения сохраненной температуры _t
    double GetStoredT() const { return _t; }
};

class Mesh    //Класс сетки
{
    std::vector<std::vector<Node*>> _mesh;
    std::vector<Node*> _hlines;
    std::vector<Node*> _vlines;
    Object& _obj;
    double _step;
    void LinkX();
    void LinkY();
    void Delnode(int, int);
    void Adapt();
    
public:
    Mesh(Object&, double);
    void ShowLinks();
    std::vector<std::vector<Node*>>& Nodes();
    std::vector<Node*>& LineX();
    std::vector<Node*>& LineY();
    Node* FindClosestNode(double x, double y) const;
    ~Mesh();
};

class System     //Класс всей пластины
{
    Object& _obj;
    Mesh _mesh;
    double _a1;
    double _a2;
    double _step;
    double _H_coeff_robin;   
    double _T_ambient_robin;
public:
   System(Object& obj, double step = 10., double a1 = 1., double a2 = 1.,
           double H_coeff_robin = 0.0, double T_ambient_robin = 0.0) // Добавлены значения по умолчанию
        : _obj(obj), _mesh(obj, step), _a1(a1), _a2(a2), _step(step),
          _H_coeff_robin(H_coeff_robin), _T_ambient_robin(T_ambient_robin) {}

    void DefineBounds(int, int, int, int);
    std::vector<std::vector<Node*>>& Nodes();
    std::vector<Node*>& LineX();
    std::vector<Node*>& LineY();
    Mesh& GetMesh();
    double step() const;
    double a1() const;
    double a2() const;

    // === ДОБАВЛЕНЫ GET-методы для параметров Робина (опционально, если нужны извне) ===
    double GetHRobin() const { return _H_coeff_robin; }
    double GetTAmbientRobin() const { return _T_ambient_robin; }
};

class Solver    //Класс вычисления решения задачи
{
    double _dt;
    std::vector<double> ThomasMethod(std::vector<std::vector<double>>&, std::vector<double>&) const;
    void SolveLine(System&, std::vector<Node*>&) const;
    std::string _name_1;
    std::string _name_2;
    
public:
    Solver(std::string name_1, std::string name_2, double dt = 1.): _dt(dt), _name_1(name_1), _name_2(name_2) {}
    void SolveExplicit(System&, double, double, double, double&) const;
    void SolveImplicit(System&, double, double, double, double&) const;
};