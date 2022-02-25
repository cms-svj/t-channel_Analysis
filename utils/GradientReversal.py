########################################################################################

#https://github.com/fungtion/DANN/blob/master/train/main.py

from torch.autograd import Function

class GradientReversalFunction(Function):
    @staticmethod
    def forward(ctx, x, lambda_):
        ctx.lambdaSVJTagger_ = lambda_
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad_output):
        dx = ctx.lambdaSVJTagger_ * grad_output.neg()
        return dx, None

########################################################################################

#https://github.com/jvanvugt/pytorch-domain-adaptation

from torch.nn import Module
from torch.autograd import Function

"""
Gradient Reversal Layer from:
Unsupervised Domain Adaptation by Backpropagation (Ganin & Lempitsky, 2015)

Forward pass is the identity function. In the backward pass,
the upstream gradients are multiplied by -lambda (i.e. gradient is reversed)
"""
class GRF(Function):
    @staticmethod
    def forward(ctx, x, lambda_):
        ctx.lambda_ = lambda_
        return x.clone()

    @staticmethod
    def backward(ctx, grads):
        lambda_ = ctx.lambda_
        lambda_ = grads.new_tensor(lambda_)
        dx = -lambda_ * grads
        return dx, None

class GradientReversal(Module):
    def __init__(self, lambda_=1):
        super(GradientReversal, self).__init__()
        self.lambda_ = lambda_

    def forward(self, x):
        return GRF.apply(x, self.lambda_)
