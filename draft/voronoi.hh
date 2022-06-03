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
    edge *last_e = NULL;
    vertex *next_v = NULL;
    vertex(point p_value, vertex *next_v_adress):p(p_value), next_v(next_v_adress){}

    float new_p_distance(point new_p)
    { return distance(new_p, p); }
    edge *in_which_triangle(point new_p)
    {
        edge *check = last_e;
        //clockwise
        while(check)
        {  
            if(point_in_triangle(p, check->stop, check->next->stop)) break;
            else check = check->previous;             
        }
        return NULL;
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
        edge *check = last_e;
        while(check)
        {
            if(check->stop == another) break;
            else check = check->previous;
        }
        return check;
    }
    void draw_voronoi()
    {
        // triangle
        edge *start_edge = last_e;
 
        while(start_edge)
        {   
            vertex *next_left, *next_right, *left, *right;
            circle this_circle, next_circle;
            right = start_edge->stop;
            left = start_edge->next->stop;
            this_circle = make_circumcircle(p, right->p, left->p);

            vertex *outside_v = find_a_vertex_connect_two_vertex(right, left, p);
            circle outside_c = make_circumcircle(right->p, left->p, outside_v->p);
            edge *oppsite = find_edge_to_another_v(left);
            oppsite->set_voronoi_edge(this_circle.center, outside_c.center); 

            next_right = left;
            next_left = start_edge->next->next->stop;
            next_circle = make_circumcircle(p, next_right->p, next_left->p);
            start_edge->next->set_voronoi_edge(this_circle.center, next_circle.center);

            start_edge->previous;
        }
    }
};

class edge
{
public:
    vertex *start, *stop;
    edge *next = NULL;
    edge *previous = NULL;
    edge *inverse = NULL;
    voronoi_edge *ve = NULL;
    edge(vertex *v1_adress, vertex *v2_adress, edge *prev)
    : start(v1_adress), stop(v2_adress), previous(prev)
    {
        prev->next = this;
        if(!prev) v1_adress->last_e = this;
        
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
        delete(ve);
        if(previous && next->previous) 
        {
            previous->next = next;
            next->previous = previous;
        }
        else if(previous)
        {
            previous->next = next;
            start->last_e = previous;
        }
        else
        {
            start->last_e->next = next;
            next->previous = NULL;
        }
    }    

    void delete_voronoi_edge()
    {
        delete(ve);
        ve = NULL;
        inverse->ve = NULL;
    }
    
    void set_voronoi_edge(point p1, point p2)
    {
        ve = new voronoi_edge(p1, p2);
        inverse->ve = ve;
    }
};

class voronoi_edge
{
public:
    point p1, p2;
    voronoi_edge(point p1_value, point p2_value):p1(p1_value), p2(p2_value){}
};

linked_list *find_vertexes_connect_two_vertex(vertex *right, vertex *left)
{
    edge *edge_right = right->last_e;
    edge *edge_left = left->last_e;
    linked_list *head = NULL;
    while (edge_left->previous)
    {
        while (edge_right->previous)
        {
            if(edge_right->stop == edge_left->stop)
            {
                head = new linked_list(edge_right->stop, head);
            }
            edge_right = edge_right->previous;
        }
        edge_right = right->last_e;
        edge_left = edge_left->previous;
    }
    return head;
};

vertex *find_a_vertex_connect_two_vertex(vertex *right, vertex *left, point p)
{
    edge *edge_right = right->last_e;
    edge *edge_left = left->last_e;
    linked_list *head = NULL;
    while (edge_left->previous)
    {
        while (edge_right->previous)
        {
            if(edge_right->stop == edge_left->stop)
            {
                head = new linked_list(edge_right->stop, head);
            }
            edge_right = edge_right->previous;
        }
        edge_right = right->last_e;
        edge_left = edge_left->previous;
    }


    return ;
}

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

bool check_in_circumcircle(vertex *v, circle circumcircle) //pending
{
    if(v)   return (distance(v->p, circumcircle.center)-circumcircle.radius)<0;
    else    return false;

}

class linked_list
{
public:
    edge *this_edge = NULL;
    vertex *this_v = NULL;
    linked_list *next_link;
    linked_list(edge *this_adress, linked_list *next_adress)
    :this_edge(this_adress), next_link(next_adress){}
    linked_list(vertex *this_adress, linked_list *next_adress)
    :this_v(this_adress), next_link(next_adress){}
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
    // delete voronoi edge
    first->delete_voronoi_edge();
    second->delete_voronoi_edge();
    edge *opposite = first->stop->find_edge_to_another_v(second->stop);
    opposite->delete_voronoi_edge();
    // split triangle into three triangle
    edge *new_e1 = new edge(newv, nearest, NULL);
    new_e1->set_inverse(new edge(nearest, newv, first, second));
    // edge *new_e1_inverse = new edge(nearest, newv, first, second);
    // new_e1_inverse->set_inverse(new_e1);
    edge *new_e2 = new edge(newv, first->stop, new_e1);
    // edge *opposite = first->stop->find_edge_to_another_v(second->stop);
    new_e2->set_inverse(new edge(first->stop, newv, opposite, opposite->next));
    edge *new_e3 = new edge(newv, second->stop, new_e2);
    opposite = second->stop->find_edge_to_another_v(nearest);
    new_e3->set_inverse(new edge(second->stop, newv, opposite, opposite->next));
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
        vertex *connect = find_vertexes_connect_two_vertex(right, left);
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
            opposite = connect->find_edge_to_another_v(insert->next->stop);
            insert->set_inverse(new edge(connect, newv, opposite, opposite->next));
            // this edge finish test. remove if from list
            linked_list *finish = head;
            head = head->next_link;
            delete(finish);
            // and add new two edges in to pendding list
            head = new linked_list(opposite->next->inverse, head);
            head = new linked_list(opposite, head);
        }
        else
        {
            // if no, move to next edge
            linked_list *finish = head;
            head = head->next_link;
            delete(finish);
        }
    }
    // create voronoi
    newv->draw_voronoi();
    return newv;
}
