#include<stdio.h>
#include<malloc.h>
define datatype int 
typedef struct lnode
{
    datatype data;
    struct lnode *next;
}LNode, *LinkList;
int Length_LinkList(LinkList L)//单链表表长
{
    LinkList p=L;
    int j=0;
    while(p->next)
    {
        p=p->next;
        j++;
    }//p指向第j个结点
    return j;
}
LinkList Get_Linklist(LinkList L,int i)//单链表按序号查找
{
    LinkList p=L;
    int j=0;
    while(p->next!=NULL&&j<i)
    {
        p=p->next;
        j++;
    }
    if(j==i)return p;
    else return NULL;
}
LinkList Locate_LinkList(LinkList L,datatypy x)//单链表按值查找
{
    LinkList p=L->next;
    while(p!=NULL&&p->data!=x)
        p=p->next;
    return p;
}
Insert_LinkList(LinkList L,int i,datatype x)//单链表插入
{
    LinkList p,s;
    p=Get_Linklist(L,i-1);
    if(p==NULL)
    {
        printf("参数i有误")；
        return 0；
    }
    else
    {   
        s=(LinkList)malloc(sizeif(LNode));
        s->data=x;
        s->next=p->next;
        p->next=s;
        return 1;
    }
}
int Del_LinkList(LinkList L,int i)//单链表删除
{
    LinkList p,s;
    p=Get_Linklist(L,i-1);
    if(p==NULL)
    {
        printf("第i-1个结点不存在")；
        return -1；
    }
    else if(p->next=NULL)
    {
        printf("第i个结点不存在")；
        return 0；
    }
    else
    {
        s=p->next;
        p->next=s->next;
        free(s);
        return 1;
    }
}