#include "header.h"
#include <iostream>

// Функции для работы с фигурами //
size_t Form::counter_ = 0;

double Form::Function(double, double) {
    return 0;
}

std::pair<double, double> Form::Deriative(double, double) {
    return {0, 0};
}

std::pair<double, double> Form::size()
{
    return {0, 0};
}

bool Form::Inhere(double, double) {
    return false;
}

std::pair<double, double> Form::missX(double)
{
    return {0, 0};
}

std::pair<double, double> Form::missY(double)
{
    return {0, 0};
}

Form::Form() {
    id_ = counter_++;
    excluded_ = false;
}

size_t Form::Get_ID() const {
    return id_;
}

bool Form::Excluded() const {
    return excluded_;
}

int Form::GetB() {return _boundtype;}

bool Form::operator==(size_t id) const {
    return id_ == id;
}

Rectangle::Rectangle(double a, double b, double h_x, double h_y, bool excluded, int btype) : a_(a), b_(b), h_x_(h_x), h_y_(h_y) {
    excluded_ = excluded;
    _boundtype = btype;
}

std::pair<double, double> Rectangle::missX(double y)
{
    return {0.5 / h_x_ + a_, -0.5 / h_x_ + a_};
}
std::pair<double, double> Rectangle::missY(double x)
{
    return {0.5 / h_y_ + b_, -0.5 / h_y_ + b_};
}

std::pair<double, double> Rectangle::size()
{
    return {1 / h_x_, 1 / h_y_};
}

double Rectangle::Function(double x, double y) {
    return std::max(h_x_ * std::abs(x - a_), h_y_ * std::abs(y - b_));
}

std::pair<double, double> Rectangle::Deriative(double x, double y) {
    return {(h_x_ / 2) * ((x - a_) / std::abs(x - a_)), (h_y_ / 2) * ((y - b_) / std::abs(y - b_))};
}

bool Rectangle::Inhere(double x, double y) {
    return Function(x, y) <= 0.5;
}

Circle::Circle(double a, double b, double h_x, double h_y, bool excluded, int btype) : a_(a), b_(b), h_x_(h_x), h_y_(h_y) {
    excluded_ = excluded;
    _boundtype = btype;
}

std::pair<double, double> Circle::missY(double x)
{
    return {std::sqrt(1 - pow((h_x_ * (x - a_)), 2))/h_y_ + b_, -std::sqrt(1 - pow((h_x_ * (x - a_)), 2))/h_y_ + b_};
}
std::pair<double, double> Circle::missX(double y)
{
    return {std::sqrt(1 - pow((h_y_ * (y - b_)), 2))/h_x_ + a_, -std::sqrt(1 - pow((h_y_ * (y - b_)), 2))/h_x_ + a_};
}

double Circle::Function(double x, double y) {
    return pow(h_x_ * (x - a_), 2) + pow(h_y_ * (y - b_), 2);
}

std::pair<double, double> Circle::Deriative(double x, double y) {
    return {2 * h_x_ * (x - a_), 2 * h_y_ * (y - b_)};
}

std::pair<double, double> Circle::size()
{
    return {1 / h_x_, 1 / h_y_};
}

bool Circle::Inhere(double x, double y) {
    return Function(x, y) <= 1;
}

Arc::Arc(double a, double b, double h_x, double h_y, bool excluded, int btype): a_(a), b_(b), h_x_(h_x), h_y_(h_y) {
    excluded_ = excluded;
    _boundtype = btype;
}

std::pair<double, double> Arc::missY(double x)
{
    return {std::sqrt(1 - pow((h_x_ * (x - a_)), 2))/h_y_ + b_, std::sqrt(1 - pow((h_x_ * (x - a_)), 2))/h_y_ + b_};
}
std::pair<double, double> Arc::missX(double y)
{
    return {std::sqrt(1 - pow((h_y_ * (y - b_)), 2))/h_x_ + a_, std::sqrt(1 - pow((h_y_ * (y - b_)), 2))/h_x_ + a_};
}

double Arc::Function(double x, double y) {
    if (x >= a_ && y >= b_) {
        return pow(h_x_ * (x - a_), 2) + pow(h_y_ * (y - b_), 2);
    }
    return -1.0;
}

std::pair<double, double> Arc::Deriative(double x, double y) {
    if (x >= a_ && y >= b_) {
        return {2 * h_x_ * (x - a_), 2 * h_y_ * (y - b_)};
    }
    if (x < a_) {
        std::cout << "x < a\n";
    }
    if (y < b_) {
        std::cout << "y < b\n";
    }
    return {-1.0, -1.0};
}

std::pair<double, double> Arc::size()
{
    return {1 / h_x_, 1 / h_y_};
}

bool Arc::Inhere(double x, double y) {
    return Function(x, y) >= 1;
}

// Функции работы с элементами сетки //
double Node::X() const {return _x;}
double Node::Y() const {return _y;}

double Node::T() const {
    if (_btype == 1)
        return 100.;

    // if (_btype == 2) {
    //     if (!_left)
    //         if (_right)
    //             return _right->T();
    //     if (!_right)
    //         if (_left)
    //             return _left->T();
    //     if (!_above)
    //         if (_below)
    //             return _below->T();
    //     if (!_below)
    //         if (_above)
    //             return _above->T();
    //     if (_right && _left) {
    //         if (_right->IsBound())
    //             return _left->T();
    //         return _right->T();
    //     }
    //     if (_above && _below) {
    //         if (_above->IsBound())
    //             return _below->T();
    //         return _above->T();
    //     }
	// }
    if (_btype == 3) {
       double T_neighbor_val;
        double delta; // Расстояние до соседа

        // Определяем соседа и расстояние до него
        if (!_left && _right) { // Узел на левой границе
            T_neighbor_val = _right->GetStoredT(); // Используем сохраненную температуру соседа
            delta = Dist(_right);
        } else if (!_right && _left) { // Узел на правой границе
            T_neighbor_val = _left->GetStoredT();
            delta = Dist(_left);
        } else if (!_above && _below) { // Узел на верхней границе (Y растет вверх)
            T_neighbor_val = _below->GetStoredT();
            delta = Dist(_below);
        } else if (!_below && _above) { // Узел на нижней границе
            T_neighbor_val = _above->GetStoredT();
            delta = Dist(_above);
        } else {
            // Этот узел помечен как тип 3, но не находится на внешней границе
            // (имеет всех соседей) или является изолированным узлом типа 3.
            // Это может произойти, если Delnode создает внутренние границы типа 3.
            // Простая формула Робина предназначена для внешних границ.
            // Возвращаем _t как запасной вариант.
            return _t;
        }

        // Формула Робина: T_s = (T_соседа + H * delta * T__среды) / (1 + H * delta)
        // где H = _node_H_robin (это h/k), T_окр_среды = _node_T_ambient_robin

        if (delta == 0) { // Избегаем деления на ноль (маловероятно)
            return T_neighbor_val;
        }

        double H_val = _node_H_robin;
        double T_ambient_val = _node_T_ambient_robin;

        double denominator = 1.0 + H_val * delta;

        // Проверка на деление на ноль или очень малое число
        if (std::abs(denominator) < 1e-9) { // 1e-9 достаточно малое число
            // Эта ситуация может возникнуть, если H_val * delta очень близко к -1.
            // Физически, H (h/k) обычно положительно.
            // Можно вернуть T_neighbor_val или T_ambient_val как приближение, или вывести ошибку.
            return T_ambient_val; // Например, если сильная конвекция, T_s стремится к T_ambient_val
        }
        return (T_neighbor_val + H_val * delta * T_ambient_val) / denominator;
    }

    return _t;
}

double Node::Dist(const Node* to) const
{
    return std::sqrt(pow(X() - to->X(), 2) + pow(Y() - to->Y(), 2));
}

Node*& Node::l() {return _left;}
Node*& Node::r() {return _right;}
Node*& Node::u() {return _above;}
Node*& Node::d() {return _below;}

void Node::LinkX(Node* l, Node* r)
{
    _left = l;
    _right = r;
}

void Node::LinkY(Node* d, Node* u)
{
    _below = d;
    _above = u;
}

void Node::SetT(double t)
{
    _t = t;
}

bool Node::IsBound() {return _btype;}

void Node::SetB(int type) {_btype = type;}

Object::Object(): _w(0), _h(0) {}

double Object::Inhere(double x, double y) 
{
    for (auto form: forms_) 
    {
        if (form->Excluded()) 
        {
            if (form->Inhere(x, y)) 
            {
                return false;
            }
        } 
        else 
        {
            if (form->Inhere(x, y)) 
            {
                return true;
            }
        }
    }
    return false;
}

void Object::Updsize()
{
    for (auto form : forms_)
    {
        _w = std::max(_w, form->size().first);
        _h = std::max(_h, form->size().second);
    }
}

bool Object::Add_Form(const std::string &name, std::map<std::string, double> &args, bool excluded, int btype) {
    if (name == "Rectangle") {
        forms_.push_back(new Rectangle(args["a"], args["b"], args["h_x"], args["h_y"], excluded, btype));
        Updsize();
        return true;
    } else if (name == "Circle") {
        forms_.push_back(new Circle(args["a"], args["b"], args["h_x"], args["h_y"], excluded, btype));
        Updsize();
        return true;
    } else if (name == "Arc") {
        forms_.push_back(new Arc(args["a"], args["b"], args["h_x"], args["h_y"], excluded, btype));
        Updsize();
        return true;
    }
    return false;
}

std::pair<double, double> Object::Fillx(double x, double y)
{
    for (auto form: forms_) 
    {
        if (form->Inhere(x, y)) 
        {
            return form->missX(y);
        } 
    }
    return {0, 0};
}

std::pair<double, double> Object::Filly(double x, double y)
{
    for (auto form: forms_) 
    {
        if (form->Inhere(x, y)) 
        {
            return form->missY(x);
        } 
    }
    return {0, 0};
}

double Object::Width() const
{
    return _w;
}

double Object::Height() const
{
    return _h;
}

bool Object::Delete_Form(size_t id) {
    return false;
}

std::vector<size_t> Object::Get_IDs() {
    std::vector<size_t> ids;
    ids.reserve(forms_.size());
    for (auto form: forms_) {
        ids.push_back(form->Get_ID());
    }
    return ids;
}

Form* Object::Who(double x, double y)
{
    for (auto form: forms_) 
    {
        if (form->Inhere(x, y)) 
        {
            return form;
        }
    }
    return forms_.back();
}

Object::~Object() 
{
    for (auto form : forms_)
        delete form;
}

// Функции для работы с сеткой //
Mesh::Mesh(Object& obj, double step): _obj(obj), _step(step)
{
    for (double y = 0; y <= _obj.Height(); y += _step)
    {
        _mesh.push_back(std::vector<Node*>());
        for (double x = 0; x <= _obj.Width(); x += _step)
        {
            _mesh.back().push_back(new Node(x, y));
        }
    }
    LinkX();
    LinkY();
    Adapt();
}

void Mesh::LinkX()
{
    for (int i = 0; i < _mesh.size(); i++)
    {
        _mesh[i][0]->LinkX(nullptr, _mesh[i][1]);
        for (int j = 1; j < _mesh[i].size() - 1; j++)
            _mesh[i][j]->LinkX(_mesh[i][j - 1], _mesh[i][j + 1]);
        _mesh[i].back()->LinkX(_mesh[i][_mesh[i].size() - 2], nullptr);
    }
    for (int i = 0; i < _mesh.size(); i++)
        _hlines.push_back(_mesh[i][0]);
}

void Mesh::LinkY()
{
    for (int j = 0; j < _mesh[0].size(); j++)
    {
        _mesh[0][j]->LinkY(nullptr, _mesh[1][j]);
        for (int i = 1; i < _mesh.size() - 1; i++)
            _mesh[i][j]->LinkY(_mesh[i - 1][j], _mesh[i + 1][j]);
        _mesh[_mesh.size() - 1][j]->LinkY(_mesh[_mesh.size() - 2][j], nullptr);
    }
    for (int i = 0; i < _mesh[0].size(); i++)
        _vlines.push_back(_mesh[0][i]);
}

void Mesh::Adapt()
{
    for (int i = 0; i < _mesh.size(); i++)
    {
        int s = _mesh[i].size();
        for (int j = 0; j < s; j++)
        {
            if (!_obj.Inhere(_mesh[i][j]->X(), _mesh[i][j]->Y()))
            {
                Delnode(i, j);
                j--;
                s--;
            }
        }
    }
}

void Mesh::ShowLinks()
{
    for (auto line : _mesh)
    {
        for (auto node : line)
        {
            if (node->d())
                std::cout << "|  ";
        }
        std::cout << '\n';
        for (auto node : line)
        {
            if (node->l())
            {
                std::cout << '-';
            }
            std::cout << 'N';
            if (node->r())
            {
                std::cout << '-';
            }
            else
            {
                std::cout << '\n';
            }   
        }
        for (auto node : line)
        {
            if (node->u())
                std::cout << "|";
            std::cout << "  ";
        }
        std::cout << '\n';
    }
}

void Mesh::Delnode(int i, int j)
{
    Node* node = _mesh[i][j];
    double bndX1 = _obj.Fillx(node->X(), node->Y()).first;
    double bndX2 = _obj.Fillx(node->X(), node->Y()).second;
    double bndY1 = _obj.Filly(node->X(), node->Y()).first;
    double bndY2 = _obj.Filly(node->X(), node->Y()).second;
    int btype = _obj.Who(node->X(), node->Y())->GetB();
    if (node->l())
    {
        if (node->l()->X() != bndX2 && node->l()->X() != bndX1)
        {
            if (bndX1 != bndX2)
            {
                Node* left = new Node(bndX2, node->Y(), btype);
                Node* right = new Node(bndX1, node->Y(), btype);
                node->l()->r() = left;
                if (node->r())
                    node->r()->l() = right;
                left->LinkX(node->l(), right);
                right->LinkX(left, node->r());
                node->l() = right;
                _mesh[i].push_back(left);
                _mesh[i].push_back(right);
            }
            else
            {
                Node* left = new Node(bndX2, node->Y(), btype);
                node->l()->r() = left;
                if (node->r())
                    node->r()->l() = left;
                left->LinkX(node->l(), node->r());
                node->l() = left;
                _mesh[i].push_back(left);
            }
        }
        else
            node->l()->r() = node->r();
    }
    if (node->r())
    {
        node->r()->l() = node->l();
    }
    if (node->d())
    {
        if (node->d()->Y() != bndY2 && node->d()->Y() != bndY1)
        {
            if (bndY2 != bndY1)
            {
                Node* down = new Node(node->X(), bndY2, btype);
                Node* up = new Node(node->X(), bndY1, btype);
                node->d()->u() = down;
                if (node->u())
                    node->u()->d() = up;
                down->LinkY(node->d(), up);
                up->LinkY(down, node->u());
                node->d() = up;
                _mesh[i].push_back(down);
                _mesh[i].push_back(up);
            }
            else
            {
                Node* down = new Node(node->X(), bndY2, btype);
                node->d()->u() = down;
                if (node->u())
                    node->u()->d() = down;
                down->LinkY(node->d(), node->u());
                node->d() = down;
                _mesh[i].push_back(down);
            }
        }
        else
            node->d()->u() = node->u();
    }
    if (node->u())
    {
        node->u()->d() = node->d();
    }
    _mesh[i].erase(_mesh[i].begin() + j);
    delete node;
}

std::vector<std::vector<Node*>>& Mesh::Nodes() {return _mesh;}
std::vector<Node*>& Mesh::LineX()  {return _hlines;}
std::vector<Node*>& Mesh::LineY()  {return _vlines;}

Node* Mesh::FindClosestNode(double x, double y) const {
    Node* closest = nullptr;
    double minDist = std::numeric_limits<double>::max();
    for (const auto& line : _mesh) {
        for (auto node : line) {
            double dist = std::sqrt(pow(node->X() - x, 2) + pow(node->Y() - y, 2));
            if (dist < minDist) {
                minDist = dist;
                closest = node;
            }
        }
    }
    return closest;
}

Mesh::~Mesh()
{
    for (auto line : _mesh)
        for (auto node : line)
            delete node;
}

// Функции для работы со всей пластиной //
void System::DefineBounds(int l_type, int t_type, int r_type, int b_type) {
    Node* cur;

    // Нижняя граница (первая горизонтальная линия)
    if (!_mesh.LineX().empty()) {
        cur = _mesh.LineX().front();
        while (cur) {
//                     if (100 < cur->X() && cur->X() < 150) {
//   cur->SetB(2);
//         } else{
            cur->SetB(b_type);
        // }
            if (b_type == 3) {
                cur->SetRobinParameters(_H_coeff_robin, _T_ambient_robin);
            }
            cur = cur->r();
        }
    }

    // Верхняя граница (последняя горизонтальная линия)
    if (!_mesh.LineX().empty()) {
        cur = _mesh.LineX().back(); // Первый узел последней горизонтальной линии
        while (cur) {
            cur->SetB(t_type);
            if (t_type == 3) {
                cur->SetRobinParameters(_H_coeff_robin, _T_ambient_robin);
            }
            cur = cur->r();
        }
    }

    // Левая граница (первая вертикальная линия)
    if (!_mesh.LineY().empty()) {
        cur = _mesh.LineY().front();
        while (cur) {
            //         if (100 < cur->Y() && cur->Y() < 150) {
            //   cur->SetB(1);
            //         } else{
            cur->SetB(l_type);
                    // }
            if (l_type == 3) {
                cur->SetRobinParameters(_H_coeff_robin, _T_ambient_robin);
            }
            cur = cur->u();
        }
    }

    // Правая граница (последняя вертикальная линия)
    if (!_mesh.LineY().empty()) {
        cur = _mesh.LineY().back(); // Первый узел последней вертикальной линии
        while (cur) { // Итерация вверх по этой линии
//         if (100 < cur->Y() && cur->Y() < 150) {
//   cur->SetB(1);
//         } else{
            cur->SetB(r_type);
        // }

            if (r_type == 3) {
                cur->SetRobinParameters(_H_coeff_robin, _T_ambient_robin);
            }
            cur = cur->u(); // Переход к следующему узлу выше
        }
    }
}

std::vector<std::vector<Node*>>& System::Nodes() {return _mesh.Nodes();}
std::vector<Node*>& System::LineX() {return _mesh.LineX();}
std::vector<Node*>& System::LineY()  {return _mesh.LineY();}
Mesh& System::GetMesh() { return _mesh; }

double System::step() const {return _step;}
double System::a1() const {return _a1;}
double System::a2() const {return _a2;}

// Функции методов решения задачи //
void Solver::SolveImplicit(System& sys, double tstop, double x_star, double y_star, double& time_at_40) const
{
    std::ofstream ExplicitOut(_name_2);
    for (double t = 0.0; t < tstop; t += _dt)
    {
        for (int i = 1; i < sys.LineX().size() - 1; i++)
        {
            std::vector<Node*> temperature;
            Node* cur = sys.LineX()[i];
            while (cur)
            {
                if (cur->r() && cur->r()->X() - cur->X() > sys.step())
                {
                    temperature.push_back(cur);
                    SolveLine(sys, temperature);
                    temperature.clear();
                    cur = cur->r();
                }
                else
                {
                    temperature.push_back(cur);
                    cur = cur->r();
                }
            }
            SolveLine(sys, temperature);
        }
        for (int i = 1; i < sys.LineY().size() - 1; i++)
        {
            std::vector<Node*> temperature;
            Node* cur = sys.LineY()[i];
            while (cur)
            {
                if (cur->u() && cur->u()->Y() - cur->Y() > sys.step())
                {
                    temperature.push_back(cur);
                    SolveLine(sys, temperature);
                    temperature.clear();
                    cur = cur->u();
                }
                else
                {
                    temperature.push_back(cur);
                    cur = cur->u();
                }
            }
            SolveLine(sys, temperature);
        }
        Node* closest = sys.GetMesh().FindClosestNode(x_star, y_star);
        double temp = closest->T();
        if (temp >= 40.0 && time_at_40 < 0) {
            time_at_40 = t;
        }
        for (auto line : sys.Nodes()) 
        {
            for (auto node : line)
                ExplicitOut << node->X() << ' ' << node->Y() << ' ' << node->T() << '\n';
        }
        ExplicitOut << "\n\n";
    }
}

// void Solver::SolveExplicit(System& sys, double tstop, double x_star, double y_star, double& time_at_40) const {
//     std::ofstream ExplicitOut(_name_1);
//     for (double t = 0.0; t < tstop; t += _dt) {
//         for (int i = 1; i < sys.LineX().size() - 1; i++) {
//             std::vector<Node*> temperature;
//             Node* cur = sys.LineX()[i];
//             while (cur) {
//                 if (cur->r() && cur->r()->X() - cur->X() > sys.step()) {
//                     temperature.push_back(cur);
//                     SolveLine(sys, temperature);
//                     temperature.clear();
//                     cur = cur->r();
//                 }
//                 else {
//                     temperature.push_back(cur);
//                     cur = cur->r();
//                 }
//             }
//             SolveLine(sys, temperature);
//         }
//         for (int i = 1; i < sys.LineY().size() - 1; i++) {
//             std::vector<Node*> temperature;
//             Node* cur = sys.LineY()[i];
//             while (cur) {
//                 if (cur->u() && cur->u()->Y() - cur->Y() > sys.step()) {
//                     temperature.push_back(cur);
//                     SolveLine(sys, temperature);
//                     temperature.clear();
//                     cur = cur->u();
//                 }
//                 else {
//                     temperature.push_back(cur);
//                     cur = cur->u();
//                 }
//             }
//             SolveLine(sys, temperature);
//         }
//         Node* closest = sys.GetMesh().FindClosestNode(x_star, y_star);
//         double temp = closest->T();
//         if (temp >= 40.0 && time_at_40 < 0) {
//             time_at_40 = t;
//         }
//         for (auto line : sys.Nodes()) {
//             for (auto node : line)
//                 ExplicitOut << node->X() << ' ' << node->Y() << ' ' << node->T() << '\n';
//         }
//         ExplicitOut << "\n\n";
//     }
// }

void Solver::SolveExplicit(System& sys, double tstop, double x_star, double y_star, double& time_at_40) const {
    std::ofstream ExplicitOut(_name_1);
    
    // Коэффициент для явной схемы (должен быть <= 0.25 для устойчивости)
    double alpha = sys.a1() / 10.0;
    double h_space = sys.step();
    double k = alpha * _dt / (h_space * h_space);
    
    for (double t = 0.0; t < tstop; t += _dt) {
        // Сохраняем старые температуры
        std::map<Node*, double> old_temps;
        for (auto& line : sys.Nodes()) {
            for (auto node : line) {
                old_temps[node] = node->T();
            }
        }
        
        // Вычисляем новые температуры по явной схеме
        for (auto& line : sys.Nodes()) {
            for (auto node : line) {
                if (!node->IsBound()) {
                    // Только для внутренних узлов вычисляем новую температуру
                    double laplacian = 0.0;
                    int neighbors_count = 0;
                    
                    if (node->l() && old_temps.count(node->l())) {
                        laplacian += old_temps[node->l()];
                        neighbors_count++;
                    }
                    if (node->r() && old_temps.count(node->r())) {
                        laplacian += old_temps[node->r()];
                        neighbors_count++;
                    }
                    if (node->u() && old_temps.count(node->u())) {
                        laplacian += old_temps[node->u()];
                        neighbors_count++;
                    }
                    if (node->d() && old_temps.count(node->d())) {
                        laplacian += old_temps[node->d()];
                        neighbors_count++;
                    }
                    
                    if (neighbors_count > 0) {
                        laplacian -= neighbors_count * old_temps[node];
                        double new_temp = old_temps[node] + k * laplacian;
                        node->SetT(new_temp);
                    }
                }
            }
        }
        
        // Проверяем температуру в целевой точке
        Node* closest = sys.GetMesh().FindClosestNode(x_star, y_star);
        double temp = closest->T();
        if (temp >= 40.0 && time_at_40 < 0) {
            time_at_40 = t;
        }
        
        // Записываем результаты
        for (auto line : sys.Nodes()) {
            for (auto node : line)
                ExplicitOut << node->X() << ' ' << node->Y() << ' ' << node->T() << '\n';
        }
        ExplicitOut << "\n\n";
    }
}

void Solver::SolveLine(System& sys, std::vector<Node*>& n) const
{
    int size = n.size() - 2;
    double mu1 = n.front()->Dist(n[1]) / sys.step();
    double mu2 = n.back()->Dist(n[n.size() - 2]) / sys.step();
    if (mu2 == 0.)
        mu2 = .1;
    double val2 = -(2 * sys.a1()) / (pow(sys.step(), 2)) - 1 / _dt;
    double val1 = sys.a1() / (pow(sys.step(), 2));
    std::vector<std::vector<double>> next(size);
    std::vector<double> right(size);
    for (int i = 0; i < next.size(); i++)
        next[i].resize(3, 0.0);
    next[0][0] = -(2 * sys.a1()) / (mu1 * pow(sys.step(), 2)) - 1 / _dt;
    next[0][1] = (2 * sys.a1()) / ((mu1 + 1) * pow(sys.step(), 2));
    next.back()[1] = (2 * sys.a1()) / ((mu2 + 1) * pow(sys.step(), 2));
    next.back()[2] = -(2 * sys.a1()) / (mu2 * pow(sys.step(), 2)) - 1 / _dt;
    for (int i = 1; i < size - 1; i++)
    {
        next[i][0] = val1;
        next[i][1] = val2;
        next[i][2] = val1;
    }
    
    for (int i = 0; i < right.size(); i++)
        right[i] = -n[i+1]->T() / _dt;
    right.front() += -(2 * sys.a1() * n.front()->T()) / (mu1 * (mu1 + 1) * pow(sys.step(), 2));
    right.back() += -(2 * sys.a1() * n.back()->T()) / (mu2 * (mu2 + 1) * pow(sys.step(), 2));
    std::vector<double> tmps = ThomasMethod(next, right);
    for (int i = 0; i < tmps.size(); i++)
        n[i + 1]->SetT(tmps[i]);
}

// void Solver::SolveLine(System& sys, std::vector<Node*>& n) const
// {
//     int size = n.size() - 2;
    
//     // Если массив слишком мал для расчета, выходим
//     if (size <= 0) return;
    
//     double mu1 = n.front()->Dist(n[1]) / sys.step();
//     double mu2 = n.back()->Dist(n[n.size() - 2]) / sys.step();
    
//     if (mu1 == 0.) mu1 = 1.0;
//     if (mu2 == 0.) mu2 = 1.0;
    
//     // КЛЮЧЕВОЕ ИЗМЕНЕНИЕ: добавляем _dt в коэффициенты
//     double val2 = -(2 * sys.a1() * _dt) / (pow(sys.step(), 2)) - 1;
//     double val1 = (sys.a1() * _dt) / (pow(sys.step(), 2));
    
//     std::vector<std::vector<double>> next(size);
//     std::vector<double> right(size);
    
//     for (int i = 0; i < size; i++)
//         next[i].resize(3, 0.0);
    
//     // Коэффициенты с учетом _dt
//     next[0][0] = -(2 * sys.a1() * _dt) / (mu1 * pow(sys.step(), 2)) - 1;
//     next[0][1] = (2 * sys.a1() * _dt) / ((mu1 + 1) * pow(sys.step(), 2));
//     next.back()[1] = (2 * sys.a1() * _dt) / ((mu2 + 1) * pow(sys.step(), 2));
//     next.back()[2] = -(2 * sys.a1() * _dt) / (mu2 * pow(sys.step(), 2)) - 1;
    
//     for (int i = 1; i < size - 1; i++)
//     {
//         next[i][0] = val1;
//         next[i][1] = val2;
//         next[i][2] = val1;
//     }
    
//     // Правая часть также должна учитывать временной шаг
//     for (int i = 0; i < size; i++)
//         right[i] = -n[i+1]->T(); // Убрано деление на _dt
    
//     right.front() += -(2 * sys.a1() * _dt * n.front()->T()) / (mu1 * (mu1 + 1) * pow(sys.step(), 2));
//     right.back() += -(2 * sys.a1() * _dt * n.back()->T()) / (mu2 * (mu2 + 1) * pow(sys.step(), 2));
    
//     std::vector<double> tmps = ThomasMethod(next, right);
//     for (int i = 0; i < tmps.size(); i++)
//         n[i + 1]->SetT(tmps[i]);
// }

std::vector<double> Solver::ThomasMethod(std::vector<std::vector<double>>& A, std::vector<double>& b) const
{
    int row = b.size() - 1;
    std::vector<double> alph(row);
    std::vector<double> bet(row);
    std::vector<double> x(b.size());
    alph[0] = -A[0][1] / A[0][0];
    bet[0] = b[0] / A[0][0];
    for (int i = 1; i < row; i++)
    {
        double a = A[i][0];
        double b1 = A[i][1];
        double c = A[i][2];
        alph[i] = -c / (a * alph[i - 1] + b1);
        bet[i] = (b[i] - a * bet[i - 1]) / (a * alph[i - 1] + b1); 
    }
    x.back() = (b.back() - A.back()[1] * bet.back()) / (A.back()[2] + A.back()[1] * alph.back());
    for (int i = row - 1; i > -1; i--)
        x[i] = alph[i] * x[i + 1] + bet[i];
    return x;
}