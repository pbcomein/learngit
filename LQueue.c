//链队，且链表带有头结点
#include<stdio.h>
#include<malloc.h>
#define datatype int 
typedef struct node
{
    datatype data;
    struct node *next;
}QNode;
typedef struct
{
    QNode *front,*rear;
}LQueue;
LQueue *Init_LQueue()
{
    LQueue *q;
    QNode *p;
    q=(LQueue*)malloc(sizeof(LQueue));
    p=(QNode*)malloc(sizeof(QNode));
    p->next=NULL;
    q->front=p;
    q->rear=p;
    return q;
}
int empty_LQueue(LQueue *q)
{
    if(q->front==q->rear)
        return 1;
    else 
        return 0;
}
void in_LQueue(LQueue *q,datatype x)
{
    QNode *s;
    s=(QNode*)malloc(sizeof(QNode));
    s->data=x;
    s->next=NULL;
    q->rear->next=s;
    q->rear=s;
}
int out_LQueue(LQueue *q, datatype *x)
{
    QNode *s;
    if(empty_LQueue(q))
    {
        printf("队空\n");
        return 0;
    }
    else
    {
        s=q->front->next;
        q->front->next=s->next;
        *x=s->data;
        free(s);
        if(q->front->next==NULL)
            q->rear=q->front;
        return 1;
    }
    
}
int main()
{
    LQueue *L;
    int n,num,m;
    int i;
    L=Init_LQueue();
    printf("初始化完成\n");
    printf("队空：%d\n",empty_LQueue(L));
    printf("请输入入队元素个数：\n");
    scanf("%d",&n);
    printf("请输入要入队的%d个元素：\n",n);
    for(i=0;i<n;i++)
    {
        scanf("%d",&num);
        in_LQueue(L,num);
    }
    out_LQueue(L,&m);
    printf("出队的元素是%d\n",m);
    return 0;
}