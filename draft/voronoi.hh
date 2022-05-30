#include<utility>
#include<cmath>
#include<algorithm>

using namespace std;
float xmin = 0, xmax = 100, ymin = 0, ymax = 100;


class point
{
public:
    float x;
    float y;
};

class vertex
{
public:
    point p;
    edge *e = NULL;
    vertex *next_v = NULL;
    vertex(point p_value, vertex *next_v_adress):p(p_value), next_v(next_v_adress){}

    float new_p_distance(point new_p)
    { return distance(new_p, p); }
    edge *in_which_triangle(point new_p)
    {
        edge *check = e;
        while(check->next->previous)
        {
            if(point_in_triangle(p, check->stop, check->next->stop))
            { return check; }
            else
            { check = check->next; }
        }
    }
    bool point_in_triangle(point new_p, vertex *v2, vertex *v3)
    {
        float d1, d2, d3;
        d1 = sign(new_p, p, v2->p);
        d2 = sign(new_p, v2->p, v3->p);
        d3 = sign(new_p, v3->p, p);
        
        bool negative, positive;
        negative = (d1<0) || (d2<0) || (d3<0);
        positive = (d1>0) || (d2>0) || (d3>0);
        return !(negative && positive);
    }
    float sign(point p1, point p2, point p3)
    {
        return (p1.x-p3.x) * (p2.y-p3.y) - (p2.x-p3.x) * (p1.y-p3.y);
    }
    edge *find_edge_to_another_v(vertex *another)
    {
        edge *check = e;
        while(check->next->previous)
        {
            if(check->stop == another) break;
            else check = check->next;
        }
        return check;
    }
};

class edge
{
public:
    vertex *start, *stop;
    edge *next = NULL;
    edge *previous = NULL;
    edge *inverse = NULL;
    bool voronoi_finish=false;
    edge(vertex *v1_adress, vertex *v2_adress, edge *prev)
    : start(v1_adress), stop(v2_adress), previous(prev)
    {
        prev->next = this;
        if(!prev) v1_adress->e = this;
        
    }
    edge(vertex *v1_adress, vertex *v2_adress, edge *prev, edge *nex)
    : start(v1_adress), stop(v2_adress), previous(prev), next(nex)
    {
        prev->next = this;
        if(nex->previous) nex->previous = this;
    }
    void set_inverse(edge *inver)
    {
        inverse = inver;
        inver->inverse = this;
    }
    void delete_myself()
    {
        if(previous) 
        {
            previous->next = next;
            next->previous = previous;
        }
        else
        {
            start->e = next;
            next->previous = NULL;
        }
    }    
};

vertex *find_a_vertex_connect_two_vertex(vertex *right, vertex *left)
{
    edge *edge_right = right->e;
    edge *edge_left = left->e;
    bool has_one = false;
    while (edge_right->stop != edge_left->stop && edge_left->next->previous && !has_one)
    {
        while (edge_right->stop != edge_left->stop && edge_right->next->previous)
        {
            edge_right = edge_right->next;
        }

        if(edge_right->stop != edge_left->stop)
        {
            edge_left = edge_left->next;
        }
        else
        {
            has_one = true;
        }
    }
    if(has_one) return edge_right->stop;
    else return NULL;
};

float distance(point p1, point p2)
{ return sqrt(pow((p1.x)-(p2.x),2) + pow(((p1.y)-(p2.y)),2)); }

class circle
{
public:
    point center;
    float radius;
};

circle make_circumcircle(point p1, point p2, point p3)
{
    // Algorithm from O'Rourke 2ed p. 189.
    float A = p2.x - p1.x,  B = p2.y - p1.y,
    C = p3.x - p1.x,  D = p3.y - p1.y,
    E = A*(p1.x+p2.x) + B*(p1.y+p2.y),
    F = C*(p1.x+p3.x) + D*(p1.y+p3.y),
    G = 2*(A*(p3.y-p2.y) - B*(p3.x-p2.x));

    circle circumcircle;
    circumcircle.center.x = (D*E-B*F)/G;
    circumcircle.center.y = (A*F-C*E)/G;
    circumcircle.radius = (distance(p1,circumcircle.center));
    return circumcircle;
}

bool check_in_circumcircle(vertex *v, circle circumcircle)
{
    if(v)   return (distance(v->p, circumcircle.center)-circumcircle.radius)<0;
    else    return false;

}

class linked_list
{
public:
    edge *this_edge = NULL;
    linked_list *next_link;
    linked_list(edge *this_adress, linked_list *next_adress)
    :this_edge(this_adress), next_link(next_adress){}
};

vertex *insert_new_point(point p, vertex *vs)
{   
    vertex *nearest = vs;
    vertex *newv = new vertex(p, vs);
    // find nearest vertex
    float min = xmax;
    while(vs->next_v)
    {
        float d = nearest->new_p_distance(p);
        if (d < min)
        {
            min = d;
            nearest = vs;
        }
        vs = vs->next_v;
    }
    // find in which trinagle near nearest vertex
    // first: right edge, second: left edge
    edge *first = nearest->in_which_triangle(p);
    edge *second = first->next;
    // split triangle into three triangle
    edge *new_e1 = new edge(newv, nearest, NULL);
    new_e1->set_inverse(new edge(nearest, newv, first, second));
    // edge *new_e1_inverse = new edge(nearest, newv, first, second);
    // new_e1_inverse->set_inverse(new_e1);
    edge *new_e2 = new edge(newv, first->stop, new_e1);
    edge *temp = first->stop->find_edge_to_another_v(second->stop);
    new_e2->set_inverse(new edge(first->stop, newv, temp, temp->next));
    edge *new_e3 = new edge(newv, second->stop, new_e2);
    temp = second->stop->find_edge_to_another_v(nearest);
    new_e3->set_inverse(new edge(second->stop, newv, temp, temp->next));
    // pending triangle to do in-circle test
    // all trianglea have same vertex -> use the edge and next edge to reresent
    linked_list *head =  new linked_list(new_e3, NULL);
    head = new linked_list(new_e2, head);
    head = new linked_list(new_e3, head);

    while(!head->next_link)
    {
        // find circumcircle
        vertex *right = head->this_edge->start;
        vertex *left = head->this_edge->stop;
        circle circumcle = make_circumcircle(newv->p, right->p, left->p);
        vertex *connect = find_a_vertex_connect_two_vertex(right, left);
        // in-circle test
        // if there is any vertex in circumcircle, 
        // the vertex connecting two points on the bottom edge must be one of them.
        if(check_in_circumcircle(connect, circumcle))
        {
            //if it is in the circumcircle, flip the invalid edge
            edge *flip = right->find_edge_to_another_v(left);
            edge *flip_inverse = flip->inverse;
            flip->delete_myself();
            flip_inverse->delete_myself();
            delete(flip);            
            delete(flip_inverse);
            // and create a new edge connect insert vertex and connect
            edge *insert = new edge(newv, connect, head->this_edge, head->this_edge->next);
            temp = connect->find_edge_to_another_v(insert->next->stop);
            insert->set_inverse(new edge(connect, newv, temp, temp->next));
            // this edge finish test. remove if from list
            linked_list *finish = head;
            head = head->next_link;
            delete(finish);
            // and add new two edges in to pendding list
            head = new linked_list(temp->next->inverse, head);
            head = new linked_list(temp, head);
        }
        else
        {
            // if no, move to next edge
            linked_list *finish = head;
            head = head->next_link;
            delete(finish);
        }
    }
    return newv;
}

