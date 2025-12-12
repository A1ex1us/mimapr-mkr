#include <iostream>
#include <fstream>
#include <unistd.h>
#include "header.h"

// Визуализация результатов при помощи GNUPLOT //
void visualize(std::ofstream &file, std::string filename, int time_end) {
    file << "set cbrange [" << 0 << ":" << 100 << "]" << std::endl;
    file << "set size ratio " << float(400) / 500 << "\nunset key\n" << "\nset palette defined (0 0 0 1, 0.25 0 1 1, 0.5 0 1 0, 0.75 1 1 0, 1 1 0 0)\n" << std::endl;
    file << "set terminal qt font 'Arial'" << std::endl;
    file << "do for [i=0:" << time_end - 1 << "]{" << std::endl;                 
    file << "plot '" << filename << "' u 1:2:3 index i w points pt 5 palette" << std::endl;
    file << "pause " << 0.00000000010 << "}" << std::endl;
    file << "pause mouse";
}

int main()
{

    double w = 500.; // Геометрические параметры
    double h = 400.;
    double S = 100.;
    double R1 = 150.; // Радиус дуги отверстия
    double R2 = 50.; // Радиус окружности
    double XS2 = 355.;
    double YS2 = 155.;

    int l = 2; // Граничные условия: левая — теплоизоляция
    int t = 3; // верхняя — конвекция
    int r = 2; // правая — теплоизоляция
    int b = 1; // нижняя — нагрев (T = 100°C)
    int r1 = 2; // внутренняя граница (ужность) — теплоизоляция
    int r2 = 2; // дуга - не применяем ГУ

    // int l = 0; // Граничные условия: левая — теплоизоляция
    // int t = 0; // верхняя — конвекция
    // int r = 0; // правая — теплоизоляция
    // int b = 0; // нижняя — нагрев (T = 100°C)
    // int r1 = 0; // внутренняя граница (окружность) — теплоизоляция
    // int r2 = 0; // дуга - не применяем ГУ
    
    double step = 10.;
    double step2 = 5.; 


    double H_coefficient = 0.1;
    double T_ambient_for_robin = 20.0; 
    
    double time_step = 0.1; //Временной шаг
    double time_end = 100.;
    
    double C = 500.; 
    

    double x_star = 200.;
    double y_star = 100.;
    
    // Переменные для хранения времени, когда T = 40°C
    double time_at_40_explicit10 = -1.0;
    double time_at_40_explicit5 = -1.0;
    double time_at_40_implicit10 = -1.0;
    double time_at_40_implicit5 = -1.0;
    
    // Вычисление параметров отверстия и дуги //
    std::map<std::string, double> base{{"a", w / 2}, {"b", h / 2}, {"h_x", 1/ w}, {"h_y", 1/ h}};
    std::map<std::string, double> circle{{"a", XS2}, {"b", YS2}, {"h_x", 1 / R2}, {"h_y", 1 / R2}};
    std::map<std::string, double> arc{{"a", w - R1}, {"b", h - R1}, {"h_x", 1 / R1}, {"h_y", 1 / R1}};
    
    Object obj;
    obj.Add_Form("Circle", circle, true, r1);
    obj.Add_Form("Arc", arc, true, r2);
    obj.Add_Form("Rectangle", base, false, 1);

    double C_10 = 2500.0;
    double C_5 = 600.0;
    double lambda_explicit = 1.0;
    double lambda_implicit = 1.5;
    double H_explicit = 0.1;
    double H_implicit = 0.15;

    System explicit10(obj, step, C_10 , 1.0, H_explicit, T_ambient_for_robin);
    explicit10.DefineBounds(l, t, r, b); // Граничные условия для ПЛАСТИНЫ

    System explicit5(obj, step2, C_5 , 1.0, H_explicit, T_ambient_for_robin);
    explicit5.DefineBounds(l, t, r, b);

    System implicit10(obj, step, C_10, 1.0, H_implicit, T_ambient_for_robin);
    implicit10.DefineBounds(l, t, r, b);

    System implicit5(obj, step2, C_5, 1.0, H_implicit, T_ambient_for_robin);
    implicit5.DefineBounds(l, t, r, b);
    
    // Сохранение результатов вычислений и отрисовка графиков //
    Solver slv10("edumma_lab_2025_rk6_64b_shliukovap_lab1_res1.dat", "edumma_lab_2025_rk6_64b_shliukovap_lab1_res2.dat", time_step);
    slv10.SolveExplicit(explicit10, time_end, x_star, y_star, time_at_40_explicit10);
    slv10.SolveImplicit(implicit10, time_end, x_star, y_star, time_at_40_implicit10);
    
    Solver slv5("edumma_lab_2025_rk6_64b_shliukovap_lab1_res3.dat", "edumma_lab_2025_rk6_64b_shliukovap_lab1_res4.dat", time_step);
    slv5.SolveExplicit(explicit5, time_end, x_star, y_star, time_at_40_explicit5);
    slv5.SolveImplicit(implicit5, time_end, x_star, y_star, time_at_40_implicit5);
    
    
    // === Вывод времени достижения 40°C в терминал ===
    std::cout << "Точка наблюдения: (" << x_star << ", " << y_star << ")\n";
    if (time_at_40_explicit10 >= 0.0)
        std::cout << "[Explicit, h=10] t* = " << time_at_40_explicit10 << " c\n";
    else
        std::cout << "[Explicit, h=10] 40°C не достигнуто за " << time_end << " c\n";

    if (time_at_40_implicit10 >= 0.0)
        std::cout << "[Implicit, h=10] t* = " << time_at_40_implicit10 << " c\n";
    else
        std::cout << "[Implicit, h=10] 40°C не достигнуто за " << time_end << " c\n";

    if (time_at_40_explicit5 >= 0.0)
        std::cout << "[Explicit, h=5]  t* = " << time_at_40_explicit5 << " c\n";
    else
        std::cout << "[Explicit, h=5]  40°C не достигнуто за " << time_end << " c\n";

    if (time_at_40_implicit5 >= 0.0)
        std::cout << "[Implicit, h=5]  t* = " << time_at_40_implicit5 << " c\n";
    else
        std::cout << "[Implicit, h=5]  40°C не достигнуто за " << time_end << " c\n";
    std::cout << std::flush;

    std::ofstream script("explicit10.plt");
    visualize(script, "edumma_lab_2025_rk6_64b_shliukovap_lab1_res1.dat", time_end);   
    script.close();
    system("gnuplot explicit10.plt");
    
    script.open("implicit10.plt");
    visualize(script, "edumma_lab_2025_rk6_64b_shliukovap_lab1_res2.dat", time_end);
    script.close();
    system("gnuplot implicit10.plt");
    
    script.open("explicit5.plt");
    visualize(script, "edumma_lab_2025_rk6_64b_shliukovap_lab1_res3.dat", time_end);
    script.close();
    system("gnuplot explicit5.plt");
    
    script.open("implicit5.plt");
    visualize(script, "edumma_lab_2025_rk6_64b_shliukovap_lab1_res4.dat", time_end);
    script.close();
    system("gnuplot implicit5.plt");
    
    return 0;
}