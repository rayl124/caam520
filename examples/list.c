#include <stdlib.h>
#include <stdio.h>

struct node_s;

struct node_s
{
  int data;
  struct node_s *next;
};

typedef struct node_s node_t;

typedef node_t* list_t;

int list_init(list_t *list)
{
  *list = NULL;
  return 0;
}

int list_free(list_t list)
{
  node_t *node;

  node = list;

  while (node) {
    node_t *next = node->next;
    free(node);
    node = next;
  }

  return 0;
}

int list_print(const list_t list)
{
  node_t *node;

  node = list;

  while (node) {
    printf("%d\n", node->data);
    node = node->next;
  }

  return 0;
}

int list_push(list_t *list, int data)
{
  if (!*list) {
    *list = malloc(sizeof(node_t));
    if (!*list) return -1;
    (*list)->data = data;
    (*list)->next = NULL;
    return 0;
  }
  else {
    node_t *node = *list;
    while (node->next) {
      node = node->next;
    }
    node->next = malloc(sizeof(node_t));
    if (!node->next) return -1;
    node->next->data = data;
    node->next->next = NULL;
  }

  return 0;
}

int main()
{
  list_t list;
  list_init(&list);

  for (int i = 0; i < 8; i++) {
    list_push(&list, i);
  }
  list_print(list);

  list_free(list);
  return 0;
}
