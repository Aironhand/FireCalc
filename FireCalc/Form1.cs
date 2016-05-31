using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace FireCalc
{
    public partial class Form1 : Form
    {
        double res;
        //масса пожарного газа в еденицу времени, соприкасающаяся с единецей поверхности стенки горного массива, кг/с
        double m;
        //удельная теплоемкость воздуха, Дж/кг*К
        double c;
        //коэффициент теплообмена между пожарными газами и поверхностью горного массива, Вт/м^2*K
        double k;
        //коэффициент теплопроводности горного массива, Вт/м^2*K
        double lambda;
        //время, с
        double[] t;
        double delta_t=10;
        //скорость пожарных газов
        double v;
        //температура пожарных газов, К
        double[] Tn;
        double delta_Tn;
        //удельная теплоемкость пожарных газов
        double Cg;
        //длина, м
        double z;
        double delta_z=1;
        //глубина
        double r;
        double delta_r=1;
        //температура горного массива, К
        double[] T;
        double delta_T;
        //коэффициент температуропроводности горного массива, м^2/с
        double alpha;
        //степень черноты горного массива
        double epsilon;
        //степень излучения, Вт/м^2*К^4
        double sigma;
        //температура горных пород на заданной глубине, К
        double Tpor;
        //постоянная температура горных пород, К
        double Tm;
        /*температура горного массива при начальных условиях
        T(P,0)=Tpor
        T(P,t)=Tm
        */
        double Tpt;
        
        public Form1()
        {
            InitializeComponent();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            ListViewItem lvi = new ListViewItem("fdsfs");
            lvi.SubItems.Add("second");
            lvi.SubItems.Add("third");
            listView1.Items.Add(lvi);
            // listView1.Columns.Add("column_z", "r");
            // lvOut.Items.Add("fds");
            //ListViewItem lvi = new ListViewItem();
            //lvi.SubItems.Add()
            //lvOut.Columns.Add("gf");
            return;
            T = new double[101];
            Tn = new double[101];
            t = new double[101];

            //рассчет delta_Tn для первой иттерации
            delta_Tn = -k / (m * Cg) * (Tn[0] - T[0]) * delta_t * delta_z / (delta_z + v * delta_t);

            //задаем T(P,0) = Tpor
            Tpt = Tpor;

            //раccчет delta_T для первой иттерации
            //первый множитель
            double first = m * c / (k * lambda) * (delta_Tn / delta_t + v * delta_Tn / delta_z);
            //второй множитель
            double second = alpha + epsilon * sigma * sigma * (Tpt * Tpt * Tpt + (Tpt * Tpt * Tn[0])
                + Tpt*Tn[0]* Tn[0]+ Tn[0]* Tn[0]* Tn[0]);
            delta_T = first * second * delta_r;

            //задаем T(P,t) = Tm
            Tpt = Tm;

            for (uint i = 1; i<=100;i++) {

                //наращиваем t,Tn и T соответсвующими delta-ми
                t[i] += delta_t;
                Tn[i] += delta_Tn;
                T[i] += delta_T;
                //рассчет delta_Tn для следующей иттерации
                delta_Tn = -k / (m * Cg) * (Tn[0] - T[0]) * delta_t * delta_z / (delta_z + v * delta_t);

                //раccчет delta_T для следующей иттерации
                //первый множитель
                first = m * c / (k * lambda) * (delta_Tn / delta_t + v * delta_Tn / delta_z);
                //второй множитель
                second = alpha + epsilon * sigma * sigma * (Tpt * Tpt * Tpt + (Tpt * Tpt * Tn[0])
                    + Tpt * Tn[0] * Tn[0] + Tn[0] * Tn[0] * Tn[0]);
                delta_T = first * second * delta_r;
            }

        }
    }
}
