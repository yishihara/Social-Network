#include <drawing.h>
#include <GL/glut.h>
#include <GL/freeglut.h>
#include <iostream>
#include <sstream>
#include <math.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

//MENU
#define STOP (100)
#define FORCE_BASED (101)
#define FUTURE_TIME (102)
#define CLUSTER_BASED (103)
#define CLUSTER_FUTURE (110)

#define CLUSTER_COLOR (120)

#define SG (130)
#define SGBC (131)
#define SGAA (132)
#define FUTURE_SG (133)

#define SEARCH (116)

using namespace boost;

SocialGraph mgraph;


float zoom = 4.0;
float xaxis = 0.0;
float yaxis = 0.0;
int width = 800, height = 600;
int flagSpring = 0;
int node=-1;
void *font = GLUT_BITMAP_TIMES_ROMAN_24;
int flagNew = 0;
float kinetic = 0;
Graph tempg;

void cvtAxis(int x1, int y1, double &x2, double &y2){
	x2 = ((double)x1/(double)width) * (2 * (double)width/zoom) + xaxis + (-width/zoom);
	y2 = ((double)(height - y1)/(double)height) * (2 * (double)height/zoom) + yaxis + (-height/zoom);
}

void drawVertices(Nodes nodes, std::vector<double> vbc){
	for(int i=0; i<nodes.size(); i++){
			glPointSize(vbc[i]);
		//glPointSize(5);
		glColor3f(nodes[i].r, nodes[i].g, nodes[i].b);
		glBegin(GL_POINTS);
		glVertex2f(nodes[i].x, nodes[i].y);
		glEnd();
	}
}

void drawNames(Nodes nodes, facebook::Name name){
	glColor3f(0,0,0);
	for(int i=0; i<nodes.size(); i++){
		glRasterPos2f(nodes[i].x, nodes[i].y);
		glutBitmapString(font, reinterpret_cast<const unsigned char*>(name[i].c_str()));
	}
}

void drawKinetic(){
	glColor3f(1.0, 0.0, 0.0);
	std::ostringstream s;
	s << "energy=" << kinetic;
	double x,y;
	cvtAxis(width-150,height-10,x,y);
	glRasterPos2f(x,y);
	glutBitmapString(font, reinterpret_cast<const unsigned char*>(s.str().c_str()));
}

void drawEdges(Nodes nodes, Graph g, MatrixScore ebc){
	EdgeIterator e, e_end;

	glColor3d(0.0, 0.0, 0.0);
	for(tie(e, e_end) = edges(g); e != e_end; e++){
		int s = source(*e, g);
		int t = target(*e, g);
		glLineWidth(ebc[s][t]);
		//glLineWidth(1);
		glBegin(GL_LINES);
		glVertex2d(nodes[s].x, nodes[s].y);
		glVertex2d(nodes[t].x, nodes[t].y);
		glEnd();
	}
}

void drawSGEdges(Nodes nodes, Edges edges, facebook::Name name){
	MatrixScore pr = mgraph.getProbability();
	
	for(int i=0; i<edges.size(); i++){
		int s = edges[i].first;
		int t = edges[i].second;

		glColor3f(1.0, 0.0, 0.0);
		glLineWidth(5.0);
		glBegin(GL_LINES);
		glVertex2d(nodes[s].x, nodes[s].y);
		glVertex2d(nodes[t].x, nodes[t].y);
		glEnd();

		glColor3f(0,0,0);
		glRasterPos2f(nodes[t].x, nodes[t].y);
		glutBitmapString(font, reinterpret_cast<const unsigned char*>(name[t].c_str()));

	}
}

void drawPREdges(Nodes nodes, facebook::Name name){
	MatrixScore pr;
	if(flagSpring == SGBC)
	pr = mgraph.getProbability();
	else if(flagSpring == SGAA)
	pr = mgraph.getProbabilityAA();

	if(node != -1){
		for(int i=0; i<pr.size(); i++){
			if(pr[node][i] >= 1.0){
				glLineWidth(5.0);
				glColor3f(1.0, 0.0, 0.0);
				glBegin(GL_LINES);
				glVertex2d(nodes[node].x, nodes[node].y);
				glVertex2d(nodes[i].x, nodes[i].y);
				glEnd();

				glColor3f(0,0,0);
				glRasterPos2f(nodes[i].x, nodes[i].y);
				glutBitmapString(font, reinterpret_cast<const unsigned char*>(name[i].c_str()));
			}
		}
	}
}

void drawText(double x, double y, std::string text, void *ft){
	glColor3f(0,0,0.0);
	cvtAxis(10,50,x,y);
	glRasterPos2f(x,y);
	glutBitmapString(font, reinterpret_cast<const unsigned char*>(text.c_str()));
}

//display graph
void drawgraph(){
	glClear(GL_COLOR_BUFFER_BIT);

	if(flagSpring == 4)
		drawEdges(mgraph.getNodes(), mgraph.getGraph(), mgraph.getEBC());
	else if(flagSpring == 5)
		drawEdges(mgraph.getNodes(), mgraph.getGraph2(), mgraph.getEBC2());
	else
		drawEdges(mgraph.getNodes(), tempg, mgraph.getEBC());

	if(flagSpring == SGBC || flagSpring == SGAA)
	drawPREdges(mgraph.getNodes(), mgraph.getName());
	if(flagSpring == FUTURE_SG)
	drawSGEdges(mgraph.getNodes(), mgraph.getSGEdges(), mgraph.getName());

	if(flagSpring == 5)
		drawVertices(mgraph.getNodes(), mgraph.getVBC2());
	else
		drawVertices(mgraph.getNodes(), mgraph.getVBC());

	if(node != -1){
		facebook::Name name = mgraph.getName();
		Nodes nodes = mgraph.getNodes();
		drawText(nodes[node].x, nodes[node].y, name[node], font);
	}

	drawKinetic();
}

//display function
void display(void){
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(xaxis + (-width / zoom), xaxis + (width / zoom), yaxis + (-height / zoom), yaxis + (height / zoom), -1.0, 1.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	drawgraph();

	glFlush();

	glutSwapBuffers();
}

//display when idle
void idle(void){
	if(flagSpring > 0){
		if(flagSpring == SG)
			kinetic = mgraph.SpringSG(node,0);
		else if(flagSpring == SGBC)
			kinetic = mgraph.SpringSG(node,1);
		else if(flagSpring == SGAA)
			kinetic = mgraph.SpringSG(node,2);
		else if(flagSpring == FUTURE_SG)
			kinetic = mgraph.SpringSG(node,3);
		else
			kinetic = mgraph.springModel(tempg);
		glutPostRedisplay();
	}
}


int findNode(int x, int y){
	Nodes nodes = mgraph.getNodes();
	for(int i=0; i<nodes.size(); i++){
		if(abs(nodes[i].x - x) <= 1 && abs(nodes[i].y - y) <= 1)
			return i;
	}
	return -1;
}

void mouse(int button, int state, int x, int y){
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){
				double x1, y1;
				cvtAxis(x, y, x1, y1);
				node = findNode(x1, y1);
				if(node != -1){
					std::cout << mgraph.getName()[node] << std::endl;
				}
			}
			break;
		default:
			break;
	}
	glutPostRedisplay();
}

void motion(int x, int y){
	if(node != -1){
		double x1, y1;
		cvtAxis(x, y, x1, y1);
		mgraph.setNode(node, x1, y1);
	}
	glutPostRedisplay();
}

void normalkeys(unsigned char key, int x, int y){
	switch(key){
		case 'i':
			zoom += zoom*.1;
			break;
		case 'o':
			zoom -= zoom*.1;
			break;
		case 'q':
			exit(0);
			break;
		case '\033':
			exit(0);
			break;
		case 'c':
			mgraph.showDetected();
			break;
		default:
			break;
	}
	glutPostRedisplay();
}

void specialkeys(int key, int x, int y){
	switch(key){
		case GLUT_KEY_UP:
			yaxis += 20/zoom;
			break;
		case GLUT_KEY_DOWN:
			yaxis -= 20/zoom;
			break;
		case GLUT_KEY_LEFT:
			xaxis -= 20/zoom;
			break;
		case GLUT_KEY_RIGHT:
			xaxis += 20/zoom;
			break;
		default:
			break;
	}
	glutPostRedisplay();
}

void item(int mode){
	switch(mode){
		case FORCE_BASED:
			flagSpring = 1;
			std::cout << "start spring model" << std::endl;
			tempg = mgraph.getGraph();
			glutIdleFunc(idle);
			break;
		case FUTURE_TIME:
			flagSpring = 1;
			std::cout << "start spring model of future time" << std::endl;
			tempg = mgraph.getGraph2();
			glutIdleFunc(idle);
			break;
		case CLUSTER_COLOR:
			std::cout << "colored clusters" << std::endl;
			mgraph.clusterColor(mgraph.getClusterGraph());
			break;
		case CLUSTER_BASED:
			flagSpring = 4;
			std::cout << "start cluster model" << std::endl;
			tempg = mgraph.getClusterGraph();
			glutIdleFunc(idle);
			break;
		case CLUSTER_FUTURE:
			flagSpring = 5;
			std::cout << "start cluster model (future)" << std::endl;
			tempg = mgraph.getClusterGraph2();
			glutIdleFunc(idle);
			break;
		case SG:
			if(node != -1){
				flagSpring = SG;
				std::cout << "drawing graph centered " << node << std::endl;
				mgraph.selectSG(node);
				tempg = mgraph.getSG();
				glutIdleFunc(idle);
			}
			else{
				std::cout << "select node" << std::endl;
			}
			break;
		case SGBC:
			if(node != -1){
			flagSpring = SGBC;
			std::cout << "drawing graph centered on " << node << " (BC)" << std::endl;
			}
			break;
		case SGAA:
			if(node != -1){
			flagSpring = SGAA;
			std::cout << "drawing graph centered on " << node << " (AA)" << std::endl;
			}
			break;
		case FUTURE_SG:
			if(node != -1){
			flagSpring = FUTURE_SG;
			std::cout << "drawing future graph centered on " << node << std::endl;
			mgraph.selectFutureSG(node);
			tempg = mgraph.getSG2();
			glutIdleFunc(idle);
			}
			break;
		case STOP:
			flagSpring = 0;
			std::cout << "stopped" << std::endl;
			glutIdleFunc(NULL);
			break;
		case SEARCH:
			std::cout << "\nwrite name: ";
			std::string temp_name;
			std::cin >> temp_name;
			facebook::Name name = mgraph.getName();
			int i;
			for(i=0; i<name.size(); i++){
				if(temp_name == name[i]){
					node = i;
					flagSpring = 10;
					mgraph.selectSG(node);
					tempg = mgraph.getSG();
					break;
				}
			}
			if(i == name.size()){
				std::cout << "did not find name" << std::endl;
				node = -1;
				flagSpring = 0;
			}
			break;
	}
	glutPostRedisplay();
}

void init(void){
	glClearColor(1.0, 1.0, 1.0, 1.0);

	glViewport(0,0,width,height);

	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
}

//main function
int mainDrawing(int argc, char *argv[]){
	//srand((unsigned int)time(NULL));
	mgraph.initialize(argc, argv);
	mgraph.test();

	//initialize glut
	glutInit(&argc, argv);

	glutInitWindowPosition(300,200);
	glutInitWindowSize(width, height);

	//display using rgba
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	//create window
	glutCreateWindow(argv[1]);
	
	init();

	//display function
	glutDisplayFunc(display);
	//glutIdleFunc(idle);

	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutKeyboardFunc(normalkeys);
	glutSpecialFunc(specialkeys);
	
	glutCreateMenu(item);
	glutAddMenuEntry("Force-based", FORCE_BASED);
	glutAddMenuEntry("Future Time", FUTURE_TIME);
	glutAddMenuEntry("Cluster", CLUSTER_BASED);
	glutAddMenuEntry("Cluster Future", CLUSTER_FUTURE);
	glutAddMenuEntry("Color Clusters", CLUSTER_COLOR);
	glutAddMenuEntry("Draw subgraph", SG);
	glutAddMenuEntry("Draw subgraph (BC)", SGBC);
	glutAddMenuEntry("Draw subgraph (AA)", SGAA);
	glutAddMenuEntry("Draw future subgraph", FUTURE_SG);
	glutAddMenuEntry("Stop", STOP);
	glutAddMenuEntry("Search", SEARCH);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	//glutPostRedisplay();

	glutMainLoop();

	return 1;
}
